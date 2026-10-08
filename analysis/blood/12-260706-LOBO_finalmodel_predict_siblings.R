#!/usr/bin/env -S Rscript --vanilla
# ---------------------------------------------------------------------------
# Train the DEFINITIVE case-control blood classifier on ALL 71 curated samples
# (the LOBO "consensus" model) and apply it to the 9 in-house sibling PAIRS.
#
# Rationale: the LOBO cross-validation reports out-of-fold performance but never
# saves a single all-data model. Here we fit that final model once - identical
# pipeline to the definitive LOBO (gene-level, curated, DRIMSeq ~Batch+Sex+Diagnosis
# feature selection, top-50, ridge) - and use it to predict the 18 discordant-sibling
# samples, which were NEVER in the case-control cohort (independent siblings_seq batches).
#
# This is a strong age control: within each sibling pair AGE IS MATCHED (mean |gap| 2.8 y),
# and both members are young. If the classifier were merely an age detector it could not tell
# proband from sibling (both young -> both look "ASD-like") and the within-pair ranking would be
# ~chance. If it ranks the ASD proband above the unaffected sibling, it is detecting an
# ASD-specific signal beyond age. NB the ABSOLUTE prob_ASD may be biased high for both members
# (the case-control model learned ASD~young); the WITHIN-PAIR ordering is the age-controlled read.
#
# Env: conda R_downstream.
# ---------------------------------------------------------------------------
DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressMessages({
  library(DRIMSeq)
  library(glmnet)
  library(pROC)
})
Sys.setenv(
  RETICULATE_PYTHON = Sys.which("python")
)
set.seed(123)
BASE <- paste0(DATA_ROOT, "/3utr/Blood/all_till_20240724/utrcalling/")
COUNTS_H5 <- paste0(
  BASE,
  "UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad"
)
PROPS_H5 <- paste0(
  BASE,
  "UTR_count_table_genelevel_bw12_above10cts_apa_named_percentages.h5ad"
)
SHEET <- paste0(DATA_ROOT, "/notebooks/blood_samples/Classifications/2025/LOPO_siblings/samplesheet_inhouse_9pairs.csv")
OUT <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/LOBO_finalmodel_predict_siblings")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TOPK <- 50L

# densify h5ad via reticulate (anndata-R builds an invalid integer dgRMatrix); index by obs_names
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
  rownames(X) <- samples
  colnames(X) <- features
  list(
    Xfs = t(X),
    gene = as.character(reticulate::py_to_r(A$var)$utr_name),
    features = features,
    samples = samples,
    obs = reticulate::py_to_r(A$obs)
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
obs <- rc$obs
rownames(obs) <- rc$samples

# --- training set = curated 71 (case-control); test set = 18 siblings (9 pairs) ---
train_ids <- rownames(obs)[obs$status_split %in% c("Train", "Test_original")]
sheet <- read.csv(SHEET, stringsAsFactors = FALSE)
sib_ids <- sheet$sample_id
stopifnot(
  all(sib_ids %in% rc$samples),
  length(intersect(train_ids, sib_ids)) == 0
)

meta_train <- data.frame(
  sample_id = train_ids,
  Diagnosis = relevel(
    factor(
      as.character(obs[train_ids, "Diagnosis"]),
      levels = c("CTRL", "ASD")
    ),
    ref = "CTRL"
  ),
  Batch = droplevels(factor(as.character(obs[train_ids, "Batch"]))),
  Sex = factor(as.character(obs[train_ids, "Sex"])),
  stringsAsFactors = FALSE
)
rownames(meta_train) <- meta_train$sample_id
cat(
  "train (curated) n=",
  nrow(meta_train),
  " ASD/CTRL:",
  sum(meta_train$Diagnosis == "ASD"),
  "/",
  sum(meta_train$Diagnosis == "CTRL"),
  " | siblings n=",
  length(sib_ids),
  " (",
  sheet$family[!duplicated(sheet$family)] |> length(),
  "pairs )\n"
)

counts_tr <- round(rc$Xfs[, train_ids, drop = FALSE])
props <- rp$Xfs
props[is.nan(props)] <- 0

# --- DRIMSeq feature selection on the full 71 (identical to definitive LOBO) ---
cdf <- data.frame(
  gene_id = gene_id,
  feature_id = feature_id,
  counts_tr,
  check.names = FALSE
)
d <- dmDSdata(counts = cdf, samples = meta_train)
d <- dmFilter(
  d,
  min_samps_gene_expr = 7,
  min_gene_expr = 10,
  min_samps_feature_expr = 3,
  min_feature_expr = 10
)
dfull <- model.matrix(~ Batch + Sex + Diagnosis, data = DRIMSeq::samples(d))
dnull <- model.matrix(~ Batch + Sex, data = DRIMSeq::samples(d))
d <- dmPrecision(d, design = dfull, BPPARAM = BiocParallel::MulticoreParam())
d <- dmFit(d, design = dfull)
d <- dmTest(d, design = dnull)
res <- DRIMSeq::results(d, level = "feature")
res <- res[order(res$pvalue), ]
res <- res[!is.na(res$pvalue), ]
top_feats <- head(res$feature_id, TOPK)
cat(
  "DRIMSeq: tested",
  nrow(res),
  "features |",
  sum(res$adj_pvalue < 0.05, na.rm = TRUE),
  "adj_p<0.05 | selected top",
  length(top_feats),
  "\n"
)

# --- fit ridge on the 71, predict the 18 siblings ---
Xtr <- as.matrix(t(props[top_feats, train_ids, drop = FALSE]))
storage.mode(Xtr) <- "double"
Xsib <- as.matrix(t(props[top_feats, sib_ids, drop = FALSE]))
storage.mode(Xsib) <- "double"
fit <- cv.glmnet(Xtr, meta_train$Diagnosis, family = "binomial", alpha = 0)
prob_sib <- as.numeric(predict(
  fit,
  newx = Xsib,
  s = "lambda.min",
  type = "response"
))

preds <- data.frame(
  family = sheet$family,
  sample_id = sib_ids,
  member = sheet$member,
  truth = sheet$Diagnosis,
  prob_ASD = prob_sib,
  stringsAsFactors = FALSE
)
preds$pred <- ifelse(preds$prob_ASD >= 0.5, "ASD", "CTRL")
preds$correct_label <- preds$pred == preds$truth

# --- within-pair ordering (age-controlled read) + pooled AUC + exact within-pair permutation ---
piv <- do.call(
  rbind,
  lapply(split(preds, preds$family), function(g) {
    a <- g$prob_ASD[g$truth == "ASD"]
    c <- g$prob_ASD[g$truth == "CTRL"]
    data.frame(family = g$family[1], ASD = a, CTRL = c, correct = a > c)
  })
)
k <- sum(piv$correct)
n <- nrow(piv)
binom_p <- binom.test(k, n, 0.5, alternative = "greater")$p.value
ro <- roc(
  preds$truth,
  preds$prob_ASD,
  levels = c("CTRL", "ASD"),
  direction = "<",
  quiet = TRUE
)
ci <- ci.auc(ro)
# exact within-pair sign-flip permutation for pooled AUC
pa <- piv$ASD
pc <- piv$CTRL
perm_auc <- apply(expand.grid(rep(list(c(0, 1)), n)), 1, function(s) {
  asd <- ifelse(s == 0, pa, pc)
  ctr <- ifelse(s == 0, pc, pa)
  as.numeric(auc(roc(
    c(rep("ASD", n), rep("CTRL", n)),
    c(asd, ctr),
    levels = c("CTRL", "ASD"),
    direction = "<",
    quiet = TRUE
  )))
})
p_perm <- mean(perm_auc >= as.numeric(auc(ro)))

write.csv(preds, file.path(OUT, "sibling_preds.csv"), row.names = FALSE)
write.csv(
  data.frame(
    feature_id = top_feats,
    rank = seq_along(top_feats),
    adj_pvalue = res$adj_pvalue[match(top_feats, res$feature_id)]
  ),
  file.path(OUT, "selected_features.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    fit = fit,
    top_feats = top_feats,
    meta_train = meta_train,
    drimseq = res
  ),
  file.path(OUT, "final_model_and_features.rds")
)
writeLines(
  c(
    paste("train_n :", nrow(meta_train)),
    paste("sib_n :", length(sib_ids)),
    paste("topk :", TOPK),
    sprintf("within-pair proband>sibling: %d/%d (binom p=%.4g)", k, n, binom_p),
    sprintf(
      "pooled AUC (proband vs sibling): %.3f (95%% CI %.3f-%.3f)",
      as.numeric(auc(ro)),
      ci[1],
      ci[3]
    ),
    sprintf("exact within-pair permutation p: %.4f", p_perm),
    paste("date :", as.character(Sys.time()))
  ),
  file.path(OUT, "run_summary.txt")
)
cat(sprintf(
  "\n=== FINAL LOBO MODEL -> SIBLING PAIRS ===\nwithin-pair proband>sibling: %d/%d (binom p=%.3g)\npooled AUC: %.3f (95%% CI %.3f-%.3f) | exact within-pair permutation p=%.3f\n",
  k,
  n,
  binom_p,
  as.numeric(auc(ro)),
  ci[1],
  ci[3],
  p_perm
))
print(piv)
cat("saved to", OUT, "\n")
