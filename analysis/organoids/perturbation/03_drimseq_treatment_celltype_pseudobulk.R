#!/usr/bin/env Rscript

# =============================================================================
# Celltype-controlled treatment-effect DRIMSeq LRT.
#
# Follow-up to 02_drimseq_treatment_main_pseudobulk.R. Instead of adjusting for the
# developmental axis with a pseudotime spline, we pseudobulk by (sample x celltype)
# and put `celltype` directly in the design, so the `induced` contrast is estimated
# within cell type, which removes the cell-type composition confound (induced organoids
# are progenitor-enriched, Li et al.). This is the cell-type-controlled sensitivity
# analysis of the pooled screen.
#
#   full  ~ celltype + super_batch + induced
#   null  ~ celltype + super_batch
#
# super_batch is retained as the confound-breaker (super_batch 137 is 100% induced,
# so the induced effect is estimated from the induced cells present in batches 206/213
# within each cell type). A raw and a celltype-stratified delta_psi are both written.
#
# Usage:
#   Rscript 03_drimseq_treatment_celltype_pseudobulk.R \
#       --input /path/merged_*.h5ad --lineage inhibitory \
#       --output_dir /path/results/inhibitory/celltype_controlled \
#       --require_targeting_guide TRUE --min_cells_per_unit 10 --ncores 16
# =============================================================================

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT") # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressPackageStartupMessages({
  library(Matrix)
})
reticulate::use_condaenv("apa-r", required = TRUE)
box::use(
  DRIMSeq[
    dmDSdata,
    dmFilter,
    dmPrecision,
    dmFit,
    dmTest,
    results,
    samples,
    counts,
    common_precision,
    genewise_precision
  ],
  dplyr[arrange, mutate, left_join, group_by, summarise, n, `%>%`],
  tibble[tibble, as_tibble]
)

# ---- args ----
argv <- commandArgs(trailingOnly = TRUE)
defaults <- list(
  input = paste0(DATA_ROOT, "/3utr/organoids/knoblich/utrcalling/merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad"),
  lineage = "inhibitory",
  output_dir = NULL,
  require_targeting_guide = "TRUE",
  control_guide = "Control2",
  balanced_batches_only = "FALSE",
  min_cells_per_unit = 10L,
  min_samps_gene_expr = 5L,
  min_gene_expr = 50L,
  min_samps_feature_expr = 5L,
  min_feature_expr = 20L,
  min_units_per_cond = 3L,
  min_cond_expr = 10L,
  seed = 42L,
  ncores = max(1L, parallel::detectCores() - 1L)
)
i <- 1L
while (i <= length(argv)) {
  key <- sub("^--", "", argv[i])
  val <- argv[i + 1L]
  if (!(key %in% names(defaults))) {
    stop("Unknown arg: ", argv[i])
  }
  cur <- defaults[[key]]
  defaults[[key]] <- if (!is.null(cur) && is.integer(cur)) {
    as.integer(val)
  } else {
    val
  }
  i <- i + 2L
}
opts <- defaults
if (is.null(opts$output_dir)) {
  stop("--output_dir required")
}
dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
bpp <- BiocParallel::MulticoreParam(workers = opts$ncores)
slug <- gsub("[^A-Za-z0-9]+", "_", tolower(opts$lineage))

# ---- helpers (as in 02_drimseq_treatment_main_pseudobulk.R) ----
pseudobulk_by_group <- function(X, groups) {
  gf <- droplevels(as.factor(groups))
  G <- Matrix::fac2sparse(gf)
  pb <- Matrix::t(G %*% X)
  colnames(pb) <- levels(gf)
  pb
}
derive_super_batch <- function(x) substr(as.character(x), 1L, 3L)
to_dgCMatrix_double <- function(X) {
  if (inherits(X, "RsparseMatrix")) {
    tmp <- as(X, "TsparseMatrix")
    return(Matrix::sparseMatrix(
      i = tmp@i + 1L,
      j = tmp@j + 1L,
      x = as.numeric(tmp@x),
      dims = dim(tmp),
      dimnames = dimnames(tmp)
    ))
  }
  if (!inherits(X, "dgCMatrix")) {
    X <- as(X, "CsparseMatrix")
  }
  if (typeof(X@x) != "double") {
    X@x <- as.numeric(X@x)
  }
  X
}
condition_balance_prefilter <- function(
  counts_df,
  samples_df,
  min_units,
  min_expr
) {
  sc <- samples_df$sample_id
  isT <- as.character(samples_df$induced) %in%
    c("TRUE", "True", "true", "T", "1")
  ctT <- as.matrix(counts_df[, sc[isT], drop = FALSE])
  ctF <- as.matrix(counts_df[, sc[!isT], drop = FALSE])
  storage.mode(ctT) <- "double"
  storage.mode(ctF) <- "double"
  gid <- counts_df$gene_id
  gT <- rowsum(ctT, gid, reorder = FALSE)
  gF <- rowsum(ctF, gid, reorder = FALSE)
  keep <- (rowSums(gT >= min_expr) >= min_units) &
    (rowSums(gF >= min_expr) >= min_units)
  kg <- rownames(gT)[keep]
  counts_df[counts_df$gene_id %in% kg, , drop = FALSE]
}
compute_delta_psi <- function(counts_d, samples_df) {
  pb_cols_raw <- samples_df$sample_id
  if (all(pb_cols_raw %in% colnames(counts_d))) {
    pb_cols <- pb_cols_raw
  } else {
    pb_cols <- make.names(pb_cols_raw) # DRIMSeq prepends 'X' to names starting with a digit
    stopifnot(all(pb_cols %in% colnames(counts_d)))
  }
  gene_id <- counts_d$gene_id
  gu <- unique(gene_id)
  cm <- as.matrix(counts_d[, pb_cols, drop = FALSE])
  storage.mode(cm) <- "double"
  isT <- as.character(samples_df$induced) %in%
    c("TRUE", "True", "true", "T", "1")
  nC <- as.numeric(samples_df$n_cells)
  ct <- as.character(samples_df$celltype)
  wT <- nC * isT
  wF <- nC * !isT
  psiT <- psiF <- numeric(length(gu))
  psiT_s <- psiF_s <- numeric(length(gu))
  for (gi in seq_along(gu)) {
    rows <- which(gene_id == gu[gi])
    K <- length(rows)
    if (K < 2L) {
      psiT[gi] <- psiF[gi] <- psiT_s[gi] <- psiF_s[gi] <- NA_real_
      next
    }
    rn <- (seq_len(K) - 1) / (K - 1)
    Y <- cm[rows, , drop = FALSE]
    tot <- colSums(Y)
    keepu <- tot > 0
    if (!any(keepu)) {
      psiT[gi] <- psiF[gi] <- psiT_s[gi] <- psiF_s[gi] <- NA_real_
      next
    }
    prop <- sweep(Y[, keepu, drop = FALSE], 2L, tot[keepu], "/")
    psi_u <- as.numeric(crossprod(rn, prop))
    iT <- isT[keepu]
    wTk <- wT[keepu]
    wFk <- wF[keepu]
    ctk <- ct[keepu]
    psiT[gi] <- if (sum(wTk) > 0) sum(psi_u * wTk) / sum(wTk) else NA_real_
    psiF[gi] <- if (sum(wFk) > 0) sum(psi_u * wFk) / sum(wFk) else NA_real_
    # celltype-stratified: within-celltype induced-uninduced, averaged over celltypes present in both
    dcs <- c()
    for (cc in unique(ctk)) {
      m <- ctk == cc
      mi <- m & iT
      mf <- m & !iT
      if (any(mi) && any(mf)) {
        pT <- sum(psi_u[mi] * wTk[mi]) / sum(wTk[mi])
        pF <- sum(psi_u[mf] * wFk[mf]) / sum(wFk[mf])
        dcs <- c(dcs, pT - pF)
      }
    }
    psiT_s[gi] <- if (length(dcs)) mean(dcs) else NA_real_ # store stratified delta directly
  }
  tibble(
    gene_id = gu,
    psi_induced_True = psiT,
    psi_induced_False = psiF,
    delta_psi_treatment_raw = psiT - psiF,
    delta_psi_treatment_ct_strat = psiT_s
  )
}

# ---- load + subset ----
cat("[load] reading h5ad ...\n")
adata <- anndata::read_h5ad(opts$input)
obs <- adata$obs
var_all <- adata$var
pt <- as.numeric(obs[[paste0("dpt_pseudotime_", opts$lineage)]])
sample_id <- as.character(obs$sample_id)
induced <- as.character(obs$induced)
celltype <- as.character(obs$celltype)
super_batch <- derive_super_batch(sample_id)
keep <- !is.na(pt)
if (isTRUE(as.logical(opts$require_targeting_guide))) {
  gRNA <- as.character(obs$gRNA)
  isI <- induced %in% c("TRUE", "True", "true", "T", "1")
  drop <- isI & (is.na(gRNA) | gRNA == opts$control_guide)
  cat("  targeting-guide cleaning drops", sum(drop & keep), "induced cells\n")
  keep <- keep & !drop
}
if (isTRUE(as.logical(opts$balanced_batches_only))) {
  ind_k <- induced %in% c("TRUE", "True", "true", "T", "1")
  both_sb <- names(which(tapply(ind_k[keep], super_batch[keep], function(z) {
    length(unique(z)) > 1L
  })))
  cat(
    "  balanced_batches_only: keeping super_batches with both conditions:",
    paste(both_sb, collapse = ", "),
    "\n"
  )
  keep <- keep & (super_batch %in% both_sb)
}
cat("  cells kept:", sum(keep), "\n")
X <- to_dgCMatrix_double(adata$X)[keep, , drop = FALSE]
sample_id <- sample_id[keep]
induced <- induced[keep]
celltype <- celltype[keep]
super_batch <- super_batch[keep]
rm(adata)
gc(verbose = FALSE)

# ---- pseudobulk by (sample x celltype) ----
groups <- paste(sample_id, gsub("[^A-Za-z0-9]+", "_", celltype), sep = "__")
pb <- as.matrix(pseudobulk_by_group(X, groups))
storage.mode(pb) <- "integer"
uk <- colnames(pb)
meta <- tibble(sample_id = uk, sample_root = sub("__[^_].*$", "", uk))
# recover per-unit sample, celltype, induced, super_batch, n_cells
split_key <- strsplit(uk, "__", fixed = TRUE)
meta$sample_root <- vapply(split_key, `[`, character(1), 1L)
meta$celltype <- vapply(
  split_key,
  function(z) paste(z[-1], collapse = "__"),
  character(1)
)
lut <- unique(data.frame(
  sample_id,
  induced,
  super_batch,
  stringsAsFactors = FALSE
))
meta <- left_join(meta, lut, by = c("sample_root" = "sample_id"))
meta$n_cells <- as.integer(table(groups)[uk])
keepu <- meta$n_cells >= opts$min_cells_per_unit
pb <- pb[, keepu, drop = FALSE]
meta <- meta[keepu, ]
meta$sample_id <- uk[keepu]
cat(
  "[pseudobulk] units:",
  nrow(meta),
  "(induced ",
  sum(meta$induced %in% c("TRUE", "True")),
  " / uninduced ",
  sum(!(meta$induced %in% c("TRUE", "True"))),
  ")\n",
  sep = ""
)
print(table(meta$super_batch, meta$induced, dnn = c("super_batch", "induced")))
cat("celltypes:\n")
print(table(meta$celltype))
write.csv(
  as.data.frame(meta),
  file.path(opts$output_dir, "pseudobulk_units.csv"),
  row.names = FALSE
)

# ---- DRIMSeq ----
counts_df <- data.frame(
  gene_id = as.character(var_all$gene_name),
  feature_id = rownames(var_all),
  pb,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
samps <- as.data.frame(meta)[, c(
  "sample_id",
  "sample_root",
  "celltype",
  "induced",
  "super_batch",
  "n_cells"
)]
stopifnot(all(colnames(pb) == samps$sample_id))
counts_df <- condition_balance_prefilter(
  counts_df,
  samps,
  opts$min_units_per_cond,
  opts$min_cond_expr
)
cat(
  "[prefilter] genes:",
  length(unique(counts_df$gene_id)),
  " features:",
  nrow(counts_df),
  "\n"
)

d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samps)
d <- DRIMSeq::dmFilter(
  d,
  min_samps_gene_expr = opts$min_samps_gene_expr,
  min_gene_expr = opts$min_gene_expr,
  min_samps_feature_expr = opts$min_samps_feature_expr,
  min_feature_expr = opts$min_feature_expr
)
counts_post <- DRIMSeq::counts(d)
cat(
  "[dmFilter] genes:",
  length(unique(counts_post$gene_id)),
  " features:",
  nrow(counts_post),
  "\n"
)

s <- DRIMSeq::samples(d)
s$celltype <- factor(as.character(s$celltype))
s$super_batch <- factor(as.character(s$super_batch))
s$induced <- factor(as.character(s$induced), levels = c("FALSE", "TRUE"))
design_full <- model.matrix(~ celltype + super_batch + induced, data = s)
design_null <- model.matrix(~ celltype + super_batch, data = s)
cat(
  "[design] full cols:",
  ncol(design_full),
  " rank:",
  qr(design_full)$rank,
  " | null cols:",
  ncol(design_null),
  " | LRT df:",
  ncol(design_full) - ncol(design_null),
  "\n"
)
stopifnot(qr(design_full)$rank == ncol(design_full))

set.seed(opts$seed)
cat("[dmPrecision] ...\n")
d <- DRIMSeq::dmPrecision(
  d,
  design = design_full,
  add_uniform = TRUE,
  BPPARAM = bpp
)
gw <- d@genewise_precision
cp <- as.numeric(d@common_precision)[1]
bad <- !is.finite(gw) | gw <= 0
if (any(bad)) {
  gw[bad] <- cp
  d@genewise_precision <- gw
  cat("  repaired", sum(bad), "precisions\n")
}
cat("[dmFit] ...\n")
d <- DRIMSeq::dmFit(d, design = design_full, BPPARAM = bpp)
cat("[dmTest] ...\n")
d_test <- DRIMSeq::dmTest(d, design = design_null, BPPARAM = bpp)
gene_res <- DRIMSeq::results(d_test, level = "gene") |>
  as_tibble() |>
  arrange(pvalue)

# ---- delta_psi (raw + celltype-stratified) ----
cat("[delta_psi] ...\n")
dpsi <- tryCatch(
  compute_delta_psi(counts_post, DRIMSeq::samples(d)),
  error = function(e) {
    cat(
      "  delta_psi FAILED:",
      conditionMessage(e),
      "\n  (writing LRT results without delta_psi)\n"
    )
    NULL
  }
)
out <- gene_res %>%
  mutate(lineage = opts$lineage, test_type = "treatment_celltype_controlled")
if (!is.null(dpsi)) {
  out <- out %>% left_join(dpsi, by = "gene_id")
}
out_path <- file.path(
  opts$output_dir,
  sprintf("drimseq_gene_%s_treatment_celltype.csv", slug)
)
write.csv(as.data.frame(out), out_path, row.names = FALSE)

sig <- sum(out$adj_pvalue < 0.05, na.rm = TRUE)
cat("\n=== RESULT ===\n")
cat("  hits padj<0.05:", sig, "\n")
sg <- out[which(out$adj_pvalue < 0.05), ]
if (!is.null(dpsi)) {
  cat(
    "  among hits, raw delta_psi:      short",
    sum(sg$delta_psi_treatment_raw < 0, na.rm = TRUE),
    "long",
    sum(sg$delta_psi_treatment_raw > 0, na.rm = TRUE),
    "\n"
  )
  cat(
    "  among hits, ct-strat delta_psi: short",
    sum(sg$delta_psi_treatment_ct_strat < 0, na.rm = TRUE),
    "long",
    sum(sg$delta_psi_treatment_ct_strat > 0, na.rm = TRUE),
    "\n"
  )
  cat(
    "  genome-wide mean ct-strat delta_psi:",
    round(mean(out$delta_psi_treatment_ct_strat, na.rm = TRUE), 4),
    "\n"
  )
}
cat("  output:", out_path, "\n")
