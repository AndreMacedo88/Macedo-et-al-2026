#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# In-house within-family LOPO ASD classifier (Fig. 4F, 4G) - common precision
# -----------------------------------------------------------------------------
# Generates the leave-one-pair-out result that figures/fig4/supp_9_lopo_within_family.ipynb
# plots.
#
# Cohort: 9 in-house ASD-proband / unaffected-sibling pairs (18 samples) sequenced
# in the same batch (smartseq v3), from the excluded sibling batches (2nd/3rd_siblings_seq)
# - independent of the discovery cohort, and batch/age/genetics/household matched within pair.
# The pairing/diagnosis is the explicit samplesheet (samplesheet_inhouse_9pairs.csv), which
# cross-checks 1:1 against the h5ad obs (family, Diagnosis) and per_pair_global_dpsi.csv.
#
# Method (as plotted in supp_9; common precision only, since genewise
# is numerically fragile here - 2/9 folds hit non-finite optimisation):
#   each fold leaves out one pair (both members); paired DRIMSeq DTU
#   (~ family + Diagnosis vs null ~ family, common precision) on the remaining 8 pairs,
#   take the top-50 features by p-value (no leakage), fit ridge cv.glmnet on those features'
#   isoform proportions, predict the held-out pair's 2 samples. Family absorbs age/sex/
#   genetic background/household, so the contrast is the within-pair ASD effect.
#   Report pooled out-of-fold ROC AUC + the within-pair proband>sibling ordering (k/9).
#
# Env: conda apa-r. Usage:
#   Rscript 09-260703-LOPO_siblings_within_family_common.R
#   (optional flags: --topk 50 --seed 123 --cores 4 --outdir <dir>)
# =============================================================================
DATA_ROOT <- Sys.getenv("M26_DATA_ROOT") # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressMessages({
  library(DRIMSeq)
  library(glmnet)
  library(pROC)
})
Sys.setenv(
  RETICULATE_PYTHON = Sys.which("python")
)

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
HERE <- paste0(DATA_ROOT, "/notebooks/blood_samples/Classifications/2025/LOPO_siblings")
BASE <- paste0(DATA_ROOT, "/3utr/Blood/all_till_20240724/utrcalling/")
COUNTS_H5 <- getarg(
  "--counts",
  paste0(BASE, "UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad")
)
PROPS_H5 <- getarg(
  "--props",
  paste0(
    BASE,
    "UTR_count_table_genelevel_bw12_above10cts_apa_named_percentages.h5ad"
  )
)
SHEET <- getarg(
  "--samplesheet",
  file.path(HERE, "samplesheet_inhouse_9pairs.csv")
)
OUTDIR <- getarg(
  "--outdir",
  paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/LOPO_siblings_common")
)
TOPK <- as.integer(getarg("--topk", "50"))
CORES <- as.integer(getarg("--cores", "4"))
SEED <- as.integer(getarg("--seed", "123"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
BPPARAM <- BiocParallel::MulticoreParam(workers = CORES)

# lenient filter suited to the small within-family paired design
FILT <- list(
  min_samps_gene_expr = 6,
  min_gene_expr = 10,
  min_samps_feature_expr = 3,
  min_feature_expr = 10
)

message("=============================================================")
message(
  "In-house within-family LOPO (common precision)  topk=",
  TOPK,
  " seed=",
  SEED
)
message("  counts: ", COUNTS_H5)
message("  props : ", PROPS_H5)
message("  sheet : ", SHEET)
message("  outdir: ", OUTDIR)
message("=============================================================")

# anndata-R builds an invalid integer dgRMatrix for the sparse count h5ads -> densify in
# Python via reticulate (convert=FALSE -> toarray). Returns features x samples with names.
read_h5ad_dense <- function(path) {
  ad_py <- reticulate::import("anndata", convert = FALSE)
  sp <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  A <- ad_py$read_h5ad(path)
  Xpy <- A$X
  if (isTRUE(reticulate::py_to_r(sp$issparse(Xpy)))) {
    Xpy <- Xpy$toarray()
  }
  Xpy <- np$asarray(Xpy, dtype = "float64")
  X <- reticulate::py_to_r(Xpy)
  samples <- as.character(reticulate::py_to_r(A$obs_names$tolist()))
  features <- as.character(reticulate::py_to_r(A$var_names$tolist()))
  var_r <- reticulate::py_to_r(A$var)
  rownames(X) <- samples
  colnames(X) <- features
  list(
    Xfs = t(X),
    gene = as.character(var_r$utr_name),
    features = features,
    samples = samples
  )
}

rc <- read_h5ad_dense(COUNTS_H5)
rp <- read_h5ad_dense(PROPS_H5)
stopifnot(
  identical(rc$features, rp$features),
  identical(rc$samples, rp$samples)
)
gene_id <- rc$gene
feature_id <- rc$features

sheet <- read.csv(SHEET, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "Diagnosis", "family") %in% colnames(sheet)))
ids <- sheet$sample_id
stopifnot(all(ids %in% rc$samples)) # all 18 present
sheet <- sheet[match(ids, sheet$sample_id), , drop = FALSE]
rownames(sheet) <- sheet$sample_id

counts_mat <- round(rc$Xfs[, ids, drop = FALSE]) # DRIMSeq expects integer counts
props_mat <- rp$Xfs[, ids, drop = FALSE]
props_mat[is.nan(props_mat)] <- 0

meta_all <- data.frame(
  sample_id = sheet$sample_id,
  Diagnosis = factor(sheet$Diagnosis, levels = c("CTRL", "ASD")),
  family = factor(sheet$family),
  stringsAsFactors = FALSE
)
rownames(meta_all) <- meta_all$sample_id
families <- levels(meta_all$family)
message(
  "n pairs (families): ",
  length(families),
  " | samples: ",
  nrow(meta_all),
  " | ASD=",
  sum(meta_all$Diagnosis == "ASD"),
  " CTRL=",
  sum(meta_all$Diagnosis == "CTRL")
)

# -------- paired DRIMSeq feature selection, common precision --------
run_drimseq_paired <- function(cmat, meta) {
  set.seed(SEED) # reproducibility only: fixes DRIMSeq's random 10%-gene common-precision
  # subset (see 11-260710b). No model/logic/leakage change; default MulticoreParam kept.
  meta <- droplevels(meta)
  df <- data.frame(
    gene_id = gene_id,
    feature_id = feature_id,
    cmat,
    check.names = FALSE
  )
  d <- dmDSdata(counts = df, samples = meta)
  d <- dmFilter(
    d,
    min_samps_gene_expr = FILT$min_samps_gene_expr,
    min_gene_expr = FILT$min_gene_expr,
    min_samps_feature_expr = FILT$min_samps_feature_expr,
    min_feature_expr = FILT$min_feature_expr
  )
  design_full <- model.matrix(~ family + Diagnosis, data = samples(d))
  design_null <- model.matrix(~family, data = samples(d))
  # Common precision: trivial genewise grid (min length 3, discarded) then overwrite the
  # genewise vector with the single common estimate -> fast + robust on the tiny paired design.
  # prec_moderation="none" is required here: the default "trended" moderation runs a loess over
  # the genewise estimates that crashes ("invalid 'x'") on this tiny 9-pair design; we discard the
  # genewise values anyway (overwritten with the common estimate below), so no moderation is needed.
  d <- dmPrecision(
    d,
    design = design_full,
    prec_moderation = "none",
    prec_grid_length = 3L,
    BPPARAM = BPPARAM
  )
  gp <- genewise_precision(d)
  gp[["genewise_precision"]] <- common_precision(d)
  genewise_precision(d) <- gp
  d <- dmFit(d, design = design_full, BPPARAM = BPPARAM)
  d <- dmTest(d, design = design_null, BPPARAM = BPPARAM)
  res <- results(d, level = "feature")
  res[order(res$pvalue), ]
}

# -------- LOPO loop (leave one pair out) --------
pred_rows <- list()
feat_rows <- list()
results_list <- list()
for (fam in families) {
  t0 <- Sys.time()
  test_ids <- meta_all$sample_id[meta_all$family == fam]
  train_ids <- meta_all$sample_id[meta_all$family != fam]
  meta_train <- meta_all[train_ids, , drop = FALSE]
  meta_test <- meta_all[test_ids, , drop = FALSE]

  sel <- tryCatch(
    run_drimseq_paired(counts_mat[, train_ids, drop = FALSE], meta_train),
    error = function(e) {
      warning("fold ", fam, " DRIMSeq failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(sel) || nrow(sel) == 0) {
    message("  [", fam, "] skipped (DRIMSeq empty)")
    next
  }
  sel <- sel[!is.na(sel$pvalue), ]
  n_sig <- sum(sel$adj_pvalue < 0.05, na.rm = TRUE)
  top_feats <- head(sel$feature_id, TOPK)

  X_train <- t(props_mat[top_feats, train_ids, drop = FALSE])
  storage.mode(X_train) <- "double"
  X_test <- t(props_mat[top_feats, test_ids, drop = FALSE])
  storage.mode(X_test) <- "double"
  set.seed(SEED) # reproducibility only: fixes cv.glmnet internal-CV fold assignment
  fit <- cv.glmnet(
    X_train,
    meta_train$Diagnosis,
    family = "binomial",
    alpha = 0
  )
  prob_asd <- as.numeric(predict(
    fit,
    newx = X_test,
    s = "lambda.min",
    type = "response"
  ))

  pred_rows[[fam]] <- data.frame(
    family = fam,
    sample_id = test_ids,
    truth = as.character(meta_test$Diagnosis),
    prob_ASD = prob_asd,
    stringsAsFactors = FALSE
  )
  feat_rows[[fam]] <- data.frame(
    feature_id = top_feats,
    fold = fam,
    n_sig = n_sig,
    rank = seq_along(top_feats),
    stringsAsFactors = FALSE
  )
  results_list[[fam]] <- list(fit = fit, n_tested = nrow(sel), n_sig = n_sig)
  message(sprintf(
    "  [%s] tested=%d sig(adj<.05)=%d selected=%d  (%.1fs)",
    fam,
    nrow(sel),
    n_sig,
    length(top_feats),
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
}
preds <- do.call(rbind, pred_rows)
feats <- do.call(rbind, feat_rows)
preds$pred <- ifelse(preds$prob_ASD >= 0.5, "ASD", "CTRL")
preds$correct <- preds$pred == preds$truth

# feature stability across folds
stab <- as.data.frame(table(feats$feature_id), stringsAsFactors = FALSE)
colnames(stab) <- c("feature_id", "n_folds_selected")
stab$prop_folds_selected <- stab$n_folds_selected / length(families)
stab <- merge(
  stab,
  unique(data.frame(feature_id = feature_id, gene_id = gene_id)),
  by = "feature_id",
  all.x = TRUE
)
stab <- stab[order(-stab$prop_folds_selected), ]

# full-data DRIMSeq for interpretation (not used in the folds -> no leakage into the metric)
message("Full-data DRIMSeq (consensus features, interpretation-only) ...")
full_res <- tryCatch(
  run_drimseq_paired(counts_mat[, meta_all$sample_id, drop = FALSE], meta_all),
  error = function(e) {
    warning(conditionMessage(e))
    NULL
  }
)
if (!is.null(full_res)) {
  full_res$rank <- seq_len(nrow(full_res))
  write.csv(
    head(full_res, 300),
    file.path(OUTDIR, "consensus_features.csv"),
    row.names = FALSE
  )
  message(
    "  full-data DTU: ",
    sum(full_res$adj_pvalue < 0.05, na.rm = TRUE),
    " adj_p<0.05 | ",
    sum(full_res$adj_pvalue < 0.1, na.rm = TRUE),
    " adj_p<0.1"
  )
}

# -------- save --------
write.csv(preds, file.path(OUTDIR, "lopo_preds.csv"), row.names = FALSE)
write.csv(feats, file.path(OUTDIR, "lopo_feat_sel.csv"), row.names = FALSE)
write.csv(
  stab,
  file.path(OUTDIR, "lopo_feature_stability.csv"),
  row.names = FALSE
)
saveRDS(results_list, file.path(OUTDIR, "lopo_results.rds"))

ro <- roc(
  preds$truth,
  preds$prob_ASD,
  levels = c("CTRL", "ASD"),
  direction = "<",
  quiet = TRUE
)
ci <- ci.auc(ro)
piv <- do.call(
  rbind,
  lapply(split(preds, preds$family), function(d) {
    a <- d$prob_ASD[d$truth == "ASD"]
    c <- d$prob_ASD[d$truth == "CTRL"]
    data.frame(family = d$family[1], ASD = a, CTRL = c, correct = a > c)
  })
)
k <- sum(piv$correct)
n <- nrow(piv)
binom_p <- binom.test(k, n, 0.5, alternative = "greater")$p.value

writeLines(
  c(
    paste("counts     :", COUNTS_H5),
    paste("props      :", PROPS_H5),
    paste("samplesheet:", SHEET),
    paste("precision  : common"),
    paste("model      : ridge (cv.glmnet alpha=0, lambda.min)"),
    paste("topk       :", TOPK),
    paste("seed       :", SEED),
    paste("n_pairs    :", n),
    paste("n_samples  :", nrow(preds)),
    paste("pooled_auc :", round(as.numeric(auc(ro)), 4)),
    paste("auc_ci95   :", paste(round(ci[c(1, 3)], 4), collapse = "-")),
    paste("within_pair:", k, "/", n, " (binom p=", signif(binom_p, 3), ")"),
    paste(
      "filter     :",
      paste(names(FILT), unlist(FILT), sep = "=", collapse = ", ")
    ),
    paste("date       :", as.character(Sys.time()))
  ),
  file.path(OUTDIR, "run_config.txt")
)

message("=============================================================")
message(sprintf(
  "DONE. pooled AUC=%.3f (95%% CI %.3f-%.3f, n=%d) | within-pair proband>sibling %d/%d (binom p=%.3g)",
  as.numeric(auc(ro)),
  ci[1],
  ci[3],
  nrow(preds),
  k,
  n,
  binom_p
))
message("Outputs in: ", OUTDIR)
