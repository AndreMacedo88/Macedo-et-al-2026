#!/usr/bin/env -S Rscript --vanilla
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
set.seed(123)
OUT <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/LOBOCV_genelevel_curated_sexadj_ridge")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

run_drimseq <- function(ct, mt, gid, fid) {
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
  dfull <- model.matrix(~ Batch + Sex + Diagnosis, data = samples(d))
  dnull <- model.matrix(~ Batch + Sex, data = samples(d))
  d <- dmPrecision(d, design = dfull, BPPARAM = BiocParallel::MulticoreParam())
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
cat(
  "curated n=",
  nrow(meta),
  " ASD/CTRL:",
  sum(meta$Diagnosis == "ASD"),
  "/",
  sum(meta$Diagnosis == "CTRL"),
  "\n"
)
cdf <- as.data.frame(t(adata$X))
cdf <- cbind(utr_name = adata$var$utr_name, cdf)
cdf <- cdf %>% select(all_of(meta$sample_id), utr_name)
pdf_ <- as.data.frame(t(adata$layers[["proportions"]]))
pdf_ <- cbind(utr_name = adata$var$utr_name, pdf_)
pdf_ <- pdf_ %>% select(all_of(meta$sample_id), utr_name)
gid <- cdf$utr_name
fid <- rownames(cdf)
batches <- levels(meta$Batch)

per <- list()
preds <- list()
fsel <- list()
for (bt in batches) {
  mtr <- meta[meta$Batch != bt, ]
  mte <- meta[meta$Batch == bt, ]
  if (nrow(mte) < 2 || length(unique(mte$Diagnosis)) < 2) {
    cat("skip", bt, "\n")
    next
  }
  mtr$Batch <- droplevels(mtr$Batch)
  sel <- tryCatch(
    run_drimseq(cdf[, mtr$sample_id, drop = FALSE], mtr, gid, fid),
    error = function(e) {
      cat("DRIMSeq fail", bt, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(sel) || nrow(sel) == 0) {
    next
  }
  feats <- head(sel$feature_id, 50)
  fsel[[bt]] <- data.frame(
    feature_id = sel$feature_id,
    adj_pvalue = sel$adj_pvalue,
    selected = sel$feature_id %in% feats,
    batch = bt
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
  per[[bt]] <- data.frame(
    batch = bt,
    n_test = nrow(mte),
    n_ctrl = sum(mte$Diagnosis == "CTRL"),
    n_case = sum(mte$Diagnosis == "ASD"),
    roc_auc = m$roc_auc,
    bal_accuracy = m$bal,
    sensitivity = m$sens,
    specificity = m$spec,
    accuracy = m$acc
  )
  preds[[bt]] <- data.frame(
    batch = bt,
    sample_id = mte$sample_id,
    truth = as.character(mte$Diagnosis),
    prob_ASD = prob
  )
  cat("done", bt, " AUC=", round(m$roc_auc, 3), "\n")
}
per_df <- bind_rows(per)
pr <- bind_rows(preds)
fs_df <- bind_rows(fsel)
saveRDS(
  list(per_batch = per_df, preds = pr, feat_sel = fs_df),
  file.path(OUT, "definitive_lobo_results.rds")
)
write.csv(per_df, file.path(OUT, "per_batch_metrics.csv"), row.names = FALSE)
# Out-of-fold predictions, read by the 10-fold comparison (04) and the Fig. 4C notebook
write.csv(pr, file.path(OUT, "lobo_preds.csv"), row.names = FALSE)
ro <- roc(
  pr$truth,
  pr$prob_ASD,
  levels = c("CTRL", "ASD"),
  direction = "<",
  quiet = TRUE
)
ci <- ci.auc(ro)
cat(sprintf(
  "\n=== DEFINITIVE (gene-level, curated, Sex-adj, ridge) ===\nMEAN per-batch AUC=%.3f | POOLED AUC=%.3f (95%% CI %.3f-%.3f) n=%d\nmean bal_acc=%.3f sens=%.3f spec=%.3f\n",
  mean(per_df$roc_auc, na.rm = TRUE),
  as.numeric(auc(ro)),
  ci[1],
  ci[3],
  nrow(pr),
  mean(per_df$bal_accuracy, na.rm = TRUE),
  mean(per_df$sensitivity, na.rm = TRUE),
  mean(per_df$specificity, na.rm = TRUE)
))
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
cat(
  "features selected in ALL",
  length(unique(per_df$batch)),
  "folds:",
  sum(fs$prop == 1),
  "\n"
)
print(head(fs, 15))

# save plots
rdf <- data.frame(fpr = 1 - ro$specificities, tpr = ro$sensitivities)
rdf <- rdf[order(rdf$fpr, rdf$tpr), ] # correct ROC point order (geom_line sorts by x only -> tied-FPR sawtooth)
ggsave(
  file.path(OUT, "pooled_roc.png"),
  ggplot(rdf, aes(fpr, tpr)) +
    geom_path(color = "#2c7fb8", linewidth = 1) +
    geom_abline(linetype = "dashed", color = "grey60") +
    annotate(
      "text",
      x = 0.62,
      y = 0.1,
      label = sprintf(
        "pooled AUC=%.3f (%.2f-%.2f)",
        as.numeric(auc(ro)),
        ci[1],
        ci[3]
      )
    ) +
    labs(
      x = "FPR",
      y = "TPR",
      title = "Definitive LOBO: gene-level, curated, Sex-adj, ridge"
    ) +
    theme_bw(),
  width = 5,
  height = 5,
  dpi = 150
)
ggsave(
  file.path(OUT, "per_batch_auc.png"),
  ggplot(
    per_df,
    aes(reorder(batch, roc_auc), roc_auc, fill = n_case / (n_case + n_ctrl))
  ) +
    geom_col() +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_text(
      aes(label = sprintf("%d/%d", n_case, n_ctrl)),
      hjust = -0.1,
      size = 3
    ) +
    coord_flip(ylim = c(0, 1)) +
    scale_fill_gradient(low = "#4dac26", high = "#d01c8b", name = "ASD frac") +
    labs(x = NULL, y = "held-out AUC", title = "Per-batch AUC (ASD/CTRL)") +
    theme_bw(),
  width = 7,
  height = 4,
  dpi = 150
)
cat("PLOTS + RESULTS saved to", OUT, "\n")
