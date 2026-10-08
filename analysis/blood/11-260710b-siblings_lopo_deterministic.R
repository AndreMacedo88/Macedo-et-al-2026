#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# Siblings within-family LOPO: DETERMINISTIC, reproducible pooled AUC (2026-07-10)
# =============================================================================
# ROOT CAUSE (from 260710-siblings_lopo_reproducible.R): the 0.79/0.815/0.827 wobble is
# NOT the classifier (over 300 cv.glmnet seeds with a FIXED feature set the pooled AUC is
# constant, sd=0) and NOT any model difference (260703 == script 04 baseline on every
# setting). It is DRIMSeq's dmPrecision estimating the common precision from a RANDOM 10%
# subset of genes; that subset was drawn from an uncontrolled RNG stream, so the per-fold
# top-50 feature ranking differed run-to-run, and with n=18 held-out samples one reordered
# case/control pair moves the pooled AUC by ~1/81 ~ 0.012.
#
# MINIMAL REPRODUCIBILITY FIX (model logic UNCHANGED; NO leakage change):
#   - set.seed(SEED) as the FIRST line inside the DRIMSeq call -> fixes WHICH 10% gene
#     subset is drawn. Per-call seeding makes feature selection depend only on
#     (SEED, fold training data), so EVERY script that selects features for a fold gets
#     the identical set, independent of any other RNG the script uses.
#   - set.seed(SEED) immediately before cv.glmnet -> fixes the internal-CV fold assignment.
#     (Default nfolds kept; ridge/alpha/lambda.min unchanged. AUC is invariant to this seed
#      anyway; we seed it only for bit-for-bit reproducibility.)
# Verified: with per-call seeding, MulticoreParam is reproducible AND identical to
# SerialParam (50/50 features) -> we keep MulticoreParam, so there is NO speed penalty.
# Everything else (DRIMSeq design/filter/prec settings, TOPK=50, leave-one-pair-out split,
# train-only feature selection) is byte-identical to 260703.
#
# ANTI-CHERRY-PICK: SEED is the project-standard 123 (not searched). We also map the AUC
# over N independent stochastic DRIMSeq realisations (random subset) to SHOW the fixed-seed
# value sits inside the natural spread; the spread is written out for a QC notebook.
# Env: conda R_downstream.
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
TOPK <- 50L; CORES <- 8L; B <- 3000L; N_STOCH <- 15L; SEED <- 123L
BPPARAM <- BiocParallel::MulticoreParam(workers = CORES)
FILT <- list(min_samps_gene_expr = 6, min_gene_expr = 10, min_samps_feature_expr = 3, min_feature_expr = 10)

read_h5ad_dense <- function(path) {
  ad_py <- reticulate::import("anndata", convert = FALSE); sp <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE); A <- ad_py$read_h5ad(path); Xpy <- A$X
  if (isTRUE(reticulate::py_to_r(sp$issparse(Xpy)))) Xpy <- Xpy$toarray()
  X <- reticulate::py_to_r(np$asarray(Xpy, dtype = "float64"))
  rownames(X) <- as.character(reticulate::py_to_r(A$obs_names$tolist()))
  colnames(X) <- as.character(reticulate::py_to_r(A$var_names$tolist()))
  list(Xfs = t(X), gene = as.character(reticulate::py_to_r(A$var)$utr_name), features = colnames(X), samples = rownames(X))
}
rc <- read_h5ad_dense(COUNTS_H5); rp <- read_h5ad_dense(PROPS_H5)
stopifnot(identical(rc$features, rp$features), identical(rc$samples, rp$samples))
gene_id <- rc$gene; feature_id <- rc$features
counts_mat <- round(rc$Xfs); props_mat <- rp$Xfs; props_mat[is.nan(props_mat)] <- 0
sheet <- read.csv(SHEET, stringsAsFactors = FALSE)
meta_all <- data.frame(sample_id = sheet$sample_id, Diagnosis = factor(sheet$Diagnosis, levels = c("CTRL","ASD")),
                       family = factor(sheet$family), stringsAsFactors = FALSE)
rownames(meta_all) <- meta_all$sample_id; families <- levels(meta_all$family); y_all <- meta_all$Diagnosis

# ---- feature selection (identical to 260703) with per-call set.seed(seed) as ONLY change ----
run_drimseq_paired <- function(cmat, meta, seed = SEED) {
  set.seed(seed)                                # <-- the reproducibility fix (fixes the 10% subset)
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
auc1 <- function(y, s) as.numeric(pROC::auc(pROC::roc(y, s, levels = c("CTRL","ASD"), direction = "<", quiet = TRUE)))
lopo_from_feats <- function(feat_by_fam) {              # classifier UNCHANGED (ridge, default nfolds), just seeded
  probs <- setNames(rep(NA_real_, nrow(meta_all)), meta_all$sample_id)
  for (fam in families) {
    feats <- feat_by_fam[[fam]]
    tr <- meta_all$sample_id[meta_all$family != fam]; te <- meta_all$sample_id[meta_all$family == fam]
    Xtr <- t(props_mat[feats, tr, drop = FALSE]); Xte <- t(props_mat[feats, te, drop = FALSE])
    storage.mode(Xtr) <- "double"; storage.mode(Xte) <- "double"; yy <- droplevels(meta_all[tr, "Diagnosis"])
    set.seed(SEED)                                      # <-- fixes cv.glmnet internal folds (nfolds default)
    fit <- cv.glmnet(Xtr, yy, family = "binomial", alpha = 0)
    probs[te] <- as.numeric(predict(fit, newx = Xte, s = "lambda.min", type = "response"))
  }
  probs
}
select_all <- function(seed) { fb <- lapply(families, function(fam) {
    tr <- meta_all$sample_id[meta_all$family != fam]
    run_drimseq_paired(counts_mat[, tr, drop = FALSE], meta_all[tr, , drop = FALSE], seed = seed)
  }); names(fb) <- families; fb }
boot_ci <- function(y, p) { set.seed(1); ii <- which(y=="ASD"); jj <- which(y=="CTRL"); bs <- numeric(B)
  for (b in seq_len(B)) { s <- c(sample(ii, replace=TRUE), sample(jj, replace=TRUE))
    bs[b] <- tryCatch(auc1(y[s], p[s]), error=function(e) NA) }; quantile(bs, c(.025,.975), na.rm=TRUE) }

## ---- determinism check: same seed -> identical fold-1 ranking on repeat ----
f1a <- run_drimseq_paired(counts_mat[, meta_all$sample_id[meta_all$family != families[1]], drop=FALSE],
                          meta_all[meta_all$sample_id[meta_all$family != families[1]], , drop=FALSE], SEED)
f1b <- run_drimseq_paired(counts_mat[, meta_all$sample_id[meta_all$family != families[1]], drop=FALSE],
                          meta_all[meta_all$sample_id[meta_all$family != families[1]], , drop=FALSE], SEED)
cat(sprintf("determinism (per-call set.seed(%d), Multicore, fold 1 identical on repeat): %s\n", SEED, identical(f1a, f1b)))

## ---- DETERMINISTIC pipeline (seed = SEED everywhere) ----
det_feats <- select_all(SEED)
prob_det <- lopo_from_feats(det_feats)
a_det <- auc1(y_all, prob_det); ci_det <- boot_ci(y_all, prob_det)
del_det <- as.numeric(pROC::ci.auc(pROC::roc(y_all, prob_det, levels=c("CTRL","ASD"), direction="<", quiet=TRUE)))
wp_det <- sum(sapply(families, function(f){ p <- prob_det[meta_all$sample_id[meta_all$family==f]]
  d <- meta_all[names(p),"Diagnosis"]; p[d=="ASD"] > p[d=="CTRL"] }))
cat(sprintf("\n== DETERMINISTIC LOPO (set.seed(%d), model unchanged) ==\n", SEED))
cat(sprintf("pooled AUC = %.3f | boot95 = [%.3f, %.3f] | DeLong95 = [%.3f, %.3f] | within-pair = %d/9\n",
            a_det, ci_det[1], ci_det[2], del_det[1], del_det[3], wp_det))
# save deterministic outputs FIRST (so a stochastic-realisation failure cannot lose them)
write.csv(data.frame(auc_deterministic = a_det, boot_lo = ci_det[1], boot_hi = ci_det[2],
                     delong_lo = del_det[1], delong_hi = del_det[3], within_pair = wp_det, seed = SEED),
          file.path(OUT, "deterministic_auc_summary.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = meta_all$sample_id, truth = as.character(y_all), prob = prob_det),
          file.path(OUT, "deterministic_lopo_preds.csv"), row.names = FALSE)

## ---- honesty: AUC over N stochastic DRIMSeq realisations (random subset) ----
# tryCatch per realisation: some random 10% subsets make DRIMSeq's optim diverge; skip those.
cat(sprintf("\nmapping AUC over %d stochastic DRIMSeq realisations...\n", N_STOCH))
stoch_auc <- rep(NA_real_, N_STOCH)
for (i in seq_len(N_STOCH)) stoch_auc[i] <- tryCatch(
  auc1(y_all, lopo_from_feats(select_all(1000 + i))), error = function(e) NA_real_)
ok <- !is.na(stoch_auc)
qs <- quantile(stoch_auc[ok], c(0,.025,.5,.975,1))
cat(sprintf("stochastic AUC (%d/%d realisations OK): mean=%.3f sd=%.3f | [min %.3f, 2.5%% %.3f, median %.3f, 97.5%% %.3f, max %.3f]\n",
            sum(ok), N_STOCH, mean(stoch_auc[ok]), sd(stoch_auc[ok]), qs[1], qs[2], qs[3], qs[4], qs[5]))
cat(sprintf("fixed-seed value %.3f is at the %.0f-th percentile of the stochastic spread\n",
            a_det, 100*mean(stoch_auc[ok] <= a_det)))
write.csv(data.frame(draw = seq_len(N_STOCH), auc = stoch_auc), file.path(OUT, "drimseq_stochastic_auc.csv"), row.names = FALSE)
cat("\nsaved ->", OUT, "\n")
