#!/usr/bin/env -S Rscript --vanilla
# ---------------------------------------------------------------------------
# Stratified 10-fold CV counterpart of the final LOBO classifier.
#
# Everything is identical to the final LOBO run
#   (03-260613-LOBOCV_genelevel_curated_sexadj_ridge.R)
#   - same data     : gene-level bw12 curated cohort (status_split Train+Test_original, n=71)
#   - same features : DRIMSeq ~Batch+Sex+Diagnosis vs null ~Batch+Sex, top-50 by pvalue,
#                     re-selected inside each training fold (no leakage), genewise precision
#   - same model    : cv.glmnet ridge (alpha=0) on the selected features' proportions
#   - same metric   : pooled out-of-fold ROC AUC + 95% CI (never the per-fold mean)
# The only change is the fold scheme: leave-one-batch-out  ->  stratified (by Diagnosis)
# random 10-fold CV. Because folds now mix batches, same-batch samples appear in both train
# and test, so this is the optimistic (batch-leaky) counterpart; the gap to the LOBO AUC
# quantifies how much of the LOBO signal is batch structure vs cross-batch signal.
#
# Repeated CV (many partitions, different seeds) checks that the result is not split-specific. Each partition
# is expensive (genewise DRIMSeq ~3-4 min/fold), so partitions are computed as independent
# processes and pooled afterwards:
#   partition mode:  Rscript this.R partition <rep_id> <seed>   -> saves one partition rds
#   aggregate mode:  Rscript this.R aggregate <n_repeats>       -> pools + plots + LOBO compare
# Repeat 1 (seed 123) is the main partition. CV_WORKERS env var caps BiocParallel workers
# per process so many partitions can run concurrently without oversubscribing the machine.
# ---------------------------------------------------------------------------
DATA_ROOT <- Sys.getenv("M26_DATA_ROOT") # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressMessages({
  library(rlang)
  library(Matrix)
})
suppressMessages(box::use(
  ggplot2[...],
  DRIMSeq[dmDSdata, dmFilter, dmPrecision, dmFit, dmTest, results, samples],
  glmnet[cv.glmnet],
  dplyr[
    filter,
    select,
    all_of,
    mutate,
    group_by,
    summarise,
    arrange,
    desc,
    pull,
    bind_rows
  ],
  tibble[tibble],
  pROC[auc, roc, ci.auc]
))
suppressMessages(library(magrittr))
# standalone Rscript (--vanilla): point reticulate at the apa-r env python (has anndata)
Sys.setenv(
  RETICULATE_PYTHON = Sys.which("python")
)

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) {
  args[1]
} else {
  stop("usage: THIS.R [partition <rep_id> <seed> | aggregate <n_repeats>]")
}
K <- 10L
BASE_SEED <- 123L
WORKERS <- as.integer(Sys.getenv("CV_WORKERS", "4"))
OUT <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/CV10fold_genelevel_curated_sexadj_ridge")
PARTDIR <- file.path(OUT, "partitions")
LOBO_PREDS <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/LOBOCV_genelevel_curated_sexadj_ridge/lobo_preds.csv")
LOBO_FS <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/LOBOCV_genelevel_curated_sexadj_ridge/feature_stability.csv")
dir.create(PARTDIR, recursive = TRUE, showWarnings = FALSE)

# --- shared helpers -------------------------------------------------------
run_drimseq <- function(ct, mt, gid, fid, use_batch = TRUE) {
  cdf <- data.frame(gene_id = gid, feature_id = fid, ct)
  rownames(cdf) <- NULL
  d <- dmDSdata(counts = cdf, samples = mt)
  d <- dmFilter(
    d,
    min_samps_gene_expr = 7,
    min_gene_expr = 10,
    min_samps_feature_expr = 3,
    min_feature_expr = 10
  )
  if (use_batch) {
    dfull <- model.matrix(~ Batch + Sex + Diagnosis, data = samples(d))
    dnull <- model.matrix(~ Batch + Sex, data = samples(d))
  } else {
    dfull <- model.matrix(~ Sex + Diagnosis, data = samples(d))
    dnull <- model.matrix(~Sex, data = samples(d))
  }
  d <- dmPrecision(
    d,
    design = dfull,
    BPPARAM = BiocParallel::MulticoreParam(workers = WORKERS)
  )
  d <- dmFit(d, design = dfull)
  d <- dmTest(d, design = dnull)
  results(d, level = "feature") %>% arrange(pvalue)
}

fmet <- function(truth, prob) {
  ro <- tryCatch(
    roc(truth, prob, levels = c("CTRL", "ASD"), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  pred <- ifelse(prob >= 0.5, "ASD", "CTRL")
  tp <- sum(pred == "ASD" & truth == "ASD")
  tn <- sum(pred == "CTRL" & truth == "CTRL")
  fp <- sum(pred == "ASD" & truth == "CTRL")
  fn <- sum(pred == "CTRL" & truth == "ASD")
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA
  list(
    roc_auc = if (!is.null(ro)) as.numeric(auc(ro)) else NA,
    bal = mean(c(sens, spec), na.rm = TRUE),
    sens = sens,
    spec = spec,
    acc = (tp + tn) / length(truth)
  )
}

# stratified k-fold: fold ids assigned within each Diagnosis class -> balanced classes per fold
make_strat_folds <- function(y, k, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  for (lev in unique(y)) {
    idx <- which(y == lev)
    idx <- sample(idx)
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

load_data <- function() {
  adata <- anndata::read_h5ad(
    paste0(DATA_ROOT, "/3utr/Blood/all_till_20240724/utrcalling/UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad")
  )
  ap <- anndata::read_h5ad(
    paste0(DATA_ROOT, "/3utr/Blood/all_till_20240724/utrcalling/UTR_count_table_genelevel_bw12_above10cts_apa_named_percentages.h5ad")
  )
  ap$X[is.nan(ap$X)] <- 0
  adata$layers[["proportions"]] <- ap$X
  if (inherits(adata$X, "sparseMatrix")) {
    adata$X <- as.matrix(adata$X)
  }
  if (inherits(adata$layers[["proportions"]], "sparseMatrix")) {
    adata$layers[["proportions"]] <- as.matrix(adata$layers[["proportions"]])
  }
  keep <- adata$obs$status_split %in% c("Train", "Test_original")
  adata <- adata[keep, ]
  adata$obs$sample_id <- gsub("-", ".", adata$obs$sample_id)
  rownames(adata$obs) <- adata$obs$sample_id
  meta <- adata$obs[, c("sample_id", "Diagnosis", "Batch", "Sex")]
  meta$sample_id <- rownames(meta)
  meta$Diagnosis <- relevel(factor(meta$Diagnosis), ref = "CTRL")
  meta$Batch <- droplevels(factor(meta$Batch))
  meta$Sex <- factor(meta$Sex)
  rownames(meta) <- meta$sample_id
  cdf <- as.data.frame(t(adata$X))
  cdf <- cbind(utr_name = adata$var$utr_name, cdf)
  cdf <- cdf %>% select(all_of(meta$sample_id), utr_name)
  pdf_ <- as.data.frame(t(adata$layers[["proportions"]]))
  pdf_ <- cbind(utr_name = adata$var$utr_name, pdf_)
  pdf_ <- pdf_ %>% select(all_of(meta$sample_id), utr_name)
  list(
    meta = meta,
    cdf = cdf,
    pdf_ = pdf_,
    gid = cdf$utr_name,
    fid = rownames(cdf)
  )
}

run_one_partition <- function(D, rep_id, seed) {
  meta <- D$meta
  cdf <- D$cdf
  pdf_ <- D$pdf_
  gid <- D$gid
  fid <- D$fid
  meta$fold <- make_strat_folds(meta$Diagnosis, K, seed)
  per <- list()
  preds <- list()
  fsel <- list()
  for (k in seq_len(K)) {
    mtr <- meta[meta$fold != k, ]
    mte <- meta[meta$fold == k, ]
    mtr$Batch <- droplevels(factor(mtr$Batch))
    sel <- tryCatch(
      run_drimseq(
        cdf[, mtr$sample_id, drop = FALSE],
        mtr,
        gid,
        fid,
        use_batch = TRUE
      ),
      error = function(e) {
        cat(
          "  DRIMSeq(+Batch) fail rep",
          rep_id,
          "fold",
          k,
          ":",
          conditionMessage(e),
          "-> retry w/o Batch\n"
        )
        NULL
      }
    )
    if (is.null(sel)) {
      sel <- tryCatch(
        run_drimseq(
          cdf[, mtr$sample_id, drop = FALSE],
          mtr,
          gid,
          fid,
          use_batch = FALSE
        ),
        error = function(e) {
          cat(
            "  DRIMSeq(-Batch) fail rep",
            rep_id,
            "fold",
            k,
            ":",
            conditionMessage(e),
            "\n"
          )
          NULL
        }
      )
    }
    if (is.null(sel) || nrow(sel) == 0) {
      next
    }
    feats <- head(sel$feature_id, 50)
    fsel[[k]] <- data.frame(
      feature_id = sel$feature_id,
      adj_pvalue = sel$adj_pvalue,
      selected = sel$feature_id %in% feats,
      fold = k
    )
    Xtr <- as.matrix(t(pdf_[feats, mtr$sample_id, drop = FALSE]))
    Xte <- as.matrix(t(pdf_[feats, mte$sample_id, drop = FALSE]))
    fit <- cv.glmnet(Xtr, mtr$Diagnosis, family = "binomial", alpha = 0)
    prob <- as.numeric(predict(
      fit,
      newx = Xte,
      s = "lambda.min",
      type = "response"
    ))
    m <- fmet(as.character(mte$Diagnosis), prob)
    per[[k]] <- data.frame(
      fold = k,
      n_test = nrow(mte),
      n_ctrl = sum(mte$Diagnosis == "CTRL"),
      n_case = sum(mte$Diagnosis == "ASD"),
      roc_auc = m$roc_auc,
      bal_accuracy = m$bal,
      sensitivity = m$sens,
      specificity = m$spec,
      accuracy = m$acc
    )
    preds[[k]] <- data.frame(
      repeat_id = rep_id,
      fold = k,
      sample_id = mte$sample_id,
      truth = as.character(mte$Diagnosis),
      prob_ASD = prob
    )
    cat(
      "  rep",
      rep_id,
      "fold",
      k,
      " n=",
      nrow(mte),
      "(",
      sum(mte$Diagnosis == "ASD"),
      "ASD/",
      sum(mte$Diagnosis == "CTRL"),
      "CTRL) AUC=",
      round(m$roc_auc, 3),
      "\n"
    )
  }
  list(
    per = bind_rows(per),
    preds = bind_rows(preds),
    fsel = bind_rows(fsel),
    rep_id = rep_id,
    seed = seed
  )
}

# --- Partition mode -------------------------------------------------------
if (MODE == "partition") {
  rep_id <- as.integer(args[2])
  seed <- as.integer(args[3])
  cat("=== PARTITION rep", rep_id, "seed", seed, "(workers", WORKERS, ") ===\n")
  D <- load_data()
  if (rep_id == 1) {
    cat(
      "curated n=",
      nrow(D$meta),
      " ASD/CTRL:",
      sum(D$meta$Diagnosis == "ASD"),
      "/",
      sum(D$meta$Diagnosis == "CTRL"),
      "\n"
    )
  }
  res <- run_one_partition(D, rep_id, seed)
  saveRDS(
    res,
    file.path(PARTDIR, sprintf("rep_%02d_seed_%d.rds", rep_id, seed))
  )
  ro <- roc(
    res$preds$truth,
    res$preds$prob_ASD,
    levels = c("CTRL", "ASD"),
    direction = "<",
    quiet = TRUE
  )
  cat(sprintf(
    "PARTITION rep %d DONE pooled AUC=%.3f n=%d\n",
    rep_id,
    as.numeric(auc(ro)),
    nrow(res$preds)
  ))

  # --- Aggregate mode -------------------------------------------------------
} else if (MODE == "aggregate") {
  n_rep <- as.integer(args[2])
  files <- sort(list.files(
    PARTDIR,
    pattern = "^rep_.*\\.rds$",
    full.names = TRUE
  ))
  stopifnot(length(files) > 0)
  all_reps <- lapply(files, readRDS)
  all_reps <- all_reps[order(sapply(all_reps, function(x) x$rep_id))]
  rep_auc <- sapply(all_reps, function(x) {
    ro <- roc(
      x$preds$truth,
      x$preds$prob_ASD,
      levels = c("CTRL", "ASD"),
      direction = "<",
      quiet = TRUE
    )
    as.numeric(auc(ro))
  })
  seeds <- sapply(all_reps, function(x) x$seed)
  reps <- sapply(all_reps, function(x) x$rep_id)

  prim <- all_reps[[which(reps == 1)]]
  per_df <- prim$per
  pr <- prim$preds
  fs_df <- prim$fsel
  ro <- roc(
    pr$truth,
    pr$prob_ASD,
    levels = c("CTRL", "ASD"),
    direction = "<",
    quiet = TRUE
  )
  ci <- ci.auc(ro)
  cv_auc <- as.numeric(auc(ro))

  lb <- read.csv(LOBO_PREDS)
  ro_lobo <- roc(
    lb$truth,
    lb$prob_ASD,
    levels = c("CTRL", "ASD"),
    direction = "<",
    quiet = TRUE
  )
  ci_lobo <- ci.auc(ro_lobo)
  lobo_auc <- as.numeric(auc(ro_lobo))

  saveRDS(
    list(
      per_fold = per_df,
      preds = pr,
      feat_sel = fs_df,
      rep_auc = rep_auc,
      seeds = seeds,
      reps = reps,
      all_reps = all_reps
    ),
    file.path(OUT, "cv10_results.rds")
  )
  write.csv(pr, file.path(OUT, "cv10_preds.csv"), row.names = FALSE)
  write.csv(per_df, file.path(OUT, "per_fold_metrics.csv"), row.names = FALSE)
  write.csv(
    data.frame(repeat_id = reps, seed = seeds, pooled_auc = rep_auc),
    file.path(OUT, "repeat_auc.csv"),
    row.names = FALSE
  )
  fs <- fs_df %>%
    group_by(feature_id) %>%
    summarise(
      n_sel = sum(selected),
      prop = mean(selected),
      med_adjp = median(adj_pvalue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(prop), med_adjp)
  write.csv(fs, file.path(OUT, "feature_stability.csv"), row.names = FALSE)

  comp <- data.frame(
    scheme = c("CV10fold_stratified", "LOBO"),
    pooled_auc = c(cv_auc, lobo_auc),
    ci_low = c(ci[1], ci_lobo[1]),
    ci_high = c(ci[3], ci_lobo[3]),
    n = c(nrow(pr), nrow(lb)),
    mean_perfold_auc = c(mean(per_df$roc_auc, na.rm = TRUE), NA),
    repeated_cv_mean = c(mean(rep_auc), NA),
    repeated_cv_sd = c(sd(rep_auc), NA),
    n_partitions = c(length(rep_auc), NA)
  )
  write.csv(comp, file.path(OUT, "comparison_summary.csv"), row.names = FALSE)

  cat(sprintf(
    "\n=== 10-FOLD CV (gene-level, curated, Sex-adj, ridge) vs LOBO ===\nCV10  headline POOLED AUC=%.3f (95%% CI %.3f-%.3f) n=%d | mean per-fold=%.3f\nLOBO  POOLED AUC=%.3f (95%% CI %.3f-%.3f) n=%d\nrepeated CV (%d partitions): mean=%.3f sd=%.3f range=[%.3f, %.3f]\n",
    cv_auc,
    ci[1],
    ci[3],
    nrow(pr),
    mean(per_df$roc_auc, na.rm = TRUE),
    lobo_auc,
    ci_lobo[1],
    ci_lobo[3],
    nrow(lb),
    length(rep_auc),
    mean(rep_auc),
    sd(rep_auc),
    min(rep_auc),
    max(rep_auc)
  ))
  cat("features selected in ALL", K, "folds:", sum(fs$prop == 1), "\n")
  print(head(fs, 15))

  rdf <- data.frame(fpr = 1 - ro$specificities, tpr = ro$sensitivities)
  rdf <- rdf[order(rdf$fpr, rdf$tpr), ] # correct ROC point order (geom_line sorts by x only -> tied-FPR sawtooth)
  ggsave(
    file.path(OUT, "pooled_roc.png"),
    ggplot(rdf, aes(fpr, tpr)) +
      geom_path(color = "#d95f0e", linewidth = 1) +
      geom_abline(linetype = "dashed", color = "grey60") +
      annotate(
        "text",
        x = 0.62,
        y = 0.1,
        label = sprintf("pooled AUC=%.3f (%.2f-%.2f)", cv_auc, ci[1], ci[3])
      ) +
      labs(
        x = "FPR",
        y = "TPR",
        title = "Stratified 10-fold CV: gene-level, curated, Sex-adj, ridge"
      ) +
      theme_bw(),
    width = 5,
    height = 5,
    dpi = 150
  )
  rdf_lobo <- data.frame(
    fpr = 1 - ro_lobo$specificities,
    tpr = ro_lobo$sensitivities
  )
  rdf_lobo <- rdf_lobo[order(rdf_lobo$fpr, rdf_lobo$tpr), ] # correct ROC point order
  ov <- rbind(
    data.frame(
      fpr = rdf$fpr,
      tpr = rdf$tpr,
      scheme = sprintf("10-fold CV (AUC %.2f)", cv_auc)
    ),
    data.frame(
      fpr = rdf_lobo$fpr,
      tpr = rdf_lobo$tpr,
      scheme = sprintf("LOBO (AUC %.2f)", lobo_auc)
    )
  )
  ggsave(
    file.path(OUT, "cv_vs_lobo_roc.png"),
    ggplot(ov, aes(fpr, tpr, color = scheme)) +
      geom_path(linewidth = 1) +
      geom_abline(linetype = "dashed", color = "grey60") +
      scale_color_manual(values = c("#d95f0e", "#2c7fb8")) +
      labs(
        x = "FPR",
        y = "TPR",
        color = NULL,
        title = "Pooled ROC: 10-fold CV vs LOBO (same data/features/model)"
      ) +
      theme_bw() +
      theme(legend.position = c(0.65, 0.15)),
    width = 5.5,
    height = 5,
    dpi = 150
  )
  ggsave(
    file.path(OUT, "per_fold_auc.png"),
    ggplot(
      per_df,
      aes(
        reorder(factor(fold), roc_auc),
        roc_auc,
        fill = n_case / (n_case + n_ctrl)
      )
    ) +
      geom_col() +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_text(
        aes(label = sprintf("%d/%d", n_case, n_ctrl)),
        hjust = -0.1,
        size = 3
      ) +
      coord_flip(ylim = c(0, 1)) +
      scale_fill_gradient(
        low = "#4dac26",
        high = "#d01c8b",
        name = "ASD frac"
      ) +
      labs(
        x = "fold",
        y = "held-out AUC",
        title = "Per-fold AUC (headline partition) - noisy; pooled is the headline"
      ) +
      theme_bw(),
    width = 7,
    height = 4,
    dpi = 150
  )
  cat("AGGREGATE done. PLOTS + RESULTS saved to", OUT, "\n")
} else {
  stop("unknown mode: ", MODE)
}
