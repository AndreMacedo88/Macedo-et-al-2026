#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# Siblings within-family LOPO: reproducibility of the pooled AUC (2026-07-10)
# =============================================================================
# The supp_9 LOPO (09-260703) and the "full LOFO" (script 04 baseline) are the same model
# (paired-DRIMSeq ~family+Diagnosis top-50 per fold + ridge cv.glmnet alpha=0,
# lambda.min; roc direction="<"). The only difference is the RNG stream feeding
# cv.glmnet's internal cross-validation fold assignment, which picks lambda.min.
# With n=18 held-out samples one reordered case/control pair moves the pooled AUC
# by ~1/81 ~ 0.012, so the value wobbles ~0.79-0.82 across seeds.
#
# This script (a) confirms the DRIMSeq feature ranking is deterministic, (b) maps the
# AUC distribution over many seeds to find the central value (and where 0.79 /
# 0.815 fall), and (c) provides deterministic estimators that remove the
# seed dependence without changing the model or choosing a seed:
#   - LOOCV internal CV (nfolds = n_train): cv.glmnet foldid = 1:n -> no RNG at all.
#   - seed-averaged ensemble: mean predicted prob over many random-fold repeats.
# Both are reported with bootstrap and DeLong 95% CIs;
# no seed is hand-picked.
# Env: conda apa-r.
# =============================================================================
DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressMessages({ library(DRIMSeq); library(glmnet); library(pROC) })
Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
BASE <- paste0(DATA_ROOT, "/3utr/Blood/all_till_20240724/utrcalling/")
COUNTS_H5 <- paste0(BASE, "UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad")
PROPS_H5  <- paste0(BASE, "UTR_count_table_genelevel_bw12_above10cts_apa_named_percentages.h5ad")
HERE <- paste0(DATA_ROOT, "/notebooks/blood_samples/Classifications/2025/LOPO_siblings")
SHEET <- file.path(HERE, "samplesheet_inhouse_9pairs.csv")
OUT <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/LOPO_siblings_regularization")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TOPK <- 50L; CORES <- 4L; N_SEED <- 300L; B <- 3000L
BPPARAM <- BiocParallel::MulticoreParam(workers = CORES)
FILT <- list(min_samps_gene_expr = 6, min_gene_expr = 10, min_samps_feature_expr = 3, min_feature_expr = 10)

read_h5ad_dense <- function(path) {
  ad_py <- reticulate::import("anndata", convert = FALSE); sp <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE); A <- ad_py$read_h5ad(path); Xpy <- A$X
  if (isTRUE(reticulate::py_to_r(sp$issparse(Xpy)))) Xpy <- Xpy$toarray()
  X <- reticulate::py_to_r(np$asarray(Xpy, dtype = "float64"))
  rownames(X) <- as.character(reticulate::py_to_r(A$obs_names$tolist()))
  colnames(X) <- as.character(reticulate::py_to_r(A$var_names$tolist()))
  list(Xfs = t(X), gene = as.character(reticulate::py_to_r(A$var)$utr_name),
       features = colnames(X), samples = rownames(X))
}
rc <- read_h5ad_dense(COUNTS_H5); rp <- read_h5ad_dense(PROPS_H5)
stopifnot(identical(rc$features, rp$features), identical(rc$samples, rp$samples))
gene_id <- rc$gene; feature_id <- rc$features
counts_mat <- round(rc$Xfs); props_mat <- rp$Xfs; props_mat[is.nan(props_mat)] <- 0
sheet <- read.csv(SHEET, stringsAsFactors = FALSE)
meta_all <- data.frame(sample_id = sheet$sample_id, Diagnosis = factor(sheet$Diagnosis, levels = c("CTRL","ASD")),
                       family = factor(sheet$family), stringsAsFactors = FALSE)
rownames(meta_all) <- meta_all$sample_id; families <- levels(meta_all$family)
y_all <- meta_all$Diagnosis

run_drimseq_paired <- function(cmat, meta) {
  meta <- droplevels(meta); df <- data.frame(gene_id = gene_id, feature_id = feature_id, cmat, check.names = FALSE)
  d <- dmDSdata(counts = df, samples = meta)
  d <- dmFilter(d, min_samps_gene_expr = FILT$min_samps_gene_expr, min_gene_expr = FILT$min_gene_expr,
                min_samps_feature_expr = FILT$min_samps_feature_expr, min_feature_expr = FILT$min_feature_expr)
  df_full <- model.matrix(~ family + Diagnosis, data = samples(d)); df_null <- model.matrix(~ family, data = samples(d))
  d <- dmPrecision(d, design = df_full, prec_moderation = "none", prec_grid_length = 3L, BPPARAM = BPPARAM)
  gp <- genewise_precision(d); gp[["genewise_precision"]] <- common_precision(d); genewise_precision(d) <- gp
  d <- dmFit(d, design = df_full, BPPARAM = BPPARAM); d <- dmTest(d, design = df_null, BPPARAM = BPPARAM)
  res <- results(d, level = "feature"); res <- res[!is.na(res$pvalue), ]; head(res[order(res$pvalue), ]$feature_id, TOPK)
}
## ---- cache per-fold DRIMSeq top-50 (deterministic) ----
cat("caching per-fold DRIMSeq top-50 ...\n")
feat_cache <- list()
for (fam in families) {
  tr <- meta_all$sample_id[meta_all$family != fam]
  feat_cache[[fam]] <- run_drimseq_paired(counts_mat[, tr, drop = FALSE], meta_all[tr, , drop = FALSE])
}
# determinism check: recompute one fold and compare
chk <- run_drimseq_paired(counts_mat[, meta_all$sample_id[meta_all$family != families[1]], drop = FALSE],
                          meta_all[meta_all$sample_id[meta_all$family != families[1]], , drop = FALSE])
cat(sprintf("DRIMSeq ranking deterministic (fold 1 identical on recompute): %s\n",
            identical(chk, feat_cache[[families[1]]])))

auc1 <- function(y, s) as.numeric(pROC::auc(pROC::roc(y, s, levels = c("CTRL","ASD"), direction = "<", quiet = TRUE)))
lopo_probs <- function(nfolds = NULL) {           # one LOPO pass; nfolds=NULL -> cv.glmnet default (10)
  probs <- setNames(rep(NA_real_, nrow(meta_all)), meta_all$sample_id)
  for (fam in families) {
    feats <- feat_cache[[fam]]
    tr <- meta_all$sample_id[meta_all$family != fam]; te <- meta_all$sample_id[meta_all$family == fam]
    Xtr <- t(props_mat[feats, tr, drop = FALSE]); Xte <- t(props_mat[feats, te, drop = FALSE])
    storage.mode(Xtr) <- "double"; storage.mode(Xte) <- "double"; yy <- droplevels(meta_all[tr, "Diagnosis"])
    nf <- if (is.null(nfolds)) 10L else nfolds
    fit <- cv.glmnet(Xtr, yy, family = "binomial", alpha = 0, nfolds = nf)
    probs[te] <- as.numeric(predict(fit, newx = Xte, s = "lambda.min", type = "response"))
  }
  probs
}
boot_ci <- function(y, p) {
  ii <- which(y == "ASD"); jj <- which(y == "CTRL"); bs <- numeric(B)
  for (b in seq_len(B)) { s <- c(sample(ii, replace = TRUE), sample(jj, replace = TRUE))
    bs[b] <- tryCatch(auc1(y[s], p[s]), error = function(e) NA) }
  quantile(bs, c(.025, .975), na.rm = TRUE)
}

# ---- (a) seed distribution (default nfolds=10, as in the original runs) ----
cat(sprintf("\nseed sweep (N=%d, nfolds=10 default) ...\n", N_SEED))
seed_auc <- numeric(N_SEED); prob_sum <- setNames(rep(0, nrow(meta_all)), meta_all$sample_id)
for (s in seq_len(N_SEED)) {
  set.seed(s); p <- lopo_probs(); seed_auc[s] <- auc1(y_all, p); prob_sum <- prob_sum + p
}
qs <- quantile(seed_auc, c(0, .025, .25, .5, .75, .975, 1))
cat(sprintf("AUC over %d seeds: mean=%.3f sd=%.3f\n", N_SEED, mean(seed_auc), sd(seed_auc)))
cat("  quantiles [min,2.5,25,50,75,97.5,max]:\n  "); cat(sprintf("%.3f ", qs)); cat("\n")
cat(sprintf("  percentile of 0.790 = %.0f%% ; of 0.815 = %.0f%%\n",
            100*mean(seed_auc <= 0.790), 100*mean(seed_auc <= 0.815)))

# ---- (b) seed-averaged ensemble (deterministic given the seed set 1..N) ----
prob_ens <- prob_sum / N_SEED
a_ens <- auc1(y_all, prob_ens); ci_ens <- boot_ci(y_all, prob_ens)
del_ens <- as.numeric(pROC::ci.auc(pROC::roc(y_all, prob_ens, levels = c("CTRL","ASD"), direction = "<", quiet = TRUE)))
cat(sprintf("\nENSEMBLE (mean prob over %d seeds): AUC=%.3f  boot95=[%.3f,%.3f]  DeLong95=[%.3f,%.3f]\n",
            N_SEED, a_ens, ci_ens[1], ci_ens[2], del_ens[1], del_ens[3]))

# ---- (c) LOOCV internal CV (deterministic; no RNG) ----
set.seed(1)  # irrelevant: nfolds=n_train -> foldid fixed
prob_loo <- lopo_probs(nfolds = 16L)   # 16 training samples per fold -> LOOCV
a_loo <- auc1(y_all, prob_loo); ci_loo <- boot_ci(y_all, prob_loo)
del_loo <- as.numeric(pROC::ci.auc(pROC::roc(y_all, prob_loo, levels = c("CTRL","ASD"), direction = "<", quiet = TRUE)))
wp_loo <- sum(sapply(families, function(f){ p <- prob_loo[meta_all$sample_id[meta_all$family==f]]
  d <- meta_all[names(p),"Diagnosis"]; p[d=="ASD"] > p[d=="CTRL"] }))
cat(sprintf("LOOCV internal CV (deterministic): AUC=%.3f  boot95=[%.3f,%.3f]  DeLong95=[%.3f,%.3f]  within-pair=%d/9\n",
            a_loo, ci_loo[1], ci_loo[2], del_loo[1], del_loo[3], wp_loo))

write.csv(data.frame(seed = seq_len(N_SEED), auc = seed_auc), file.path(OUT, "seed_auc_distribution.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = meta_all$sample_id, truth = as.character(y_all),
                     prob_ensemble = prob_ens, prob_loocv = prob_loo),
          file.path(OUT, "reproducible_lopo_preds.csv"), row.names = FALSE)
cat("\nsaved -> ", OUT, "\n")
