#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# 01-261007-DTU_curated_genelevel_unadj.R
#
# Purpose
#   Gene-level DRIMSeq differential transcript (3'UTR isoform) usage, ASD vs CTRL,
#   on the curated blood cohort (71 samples: obs status_split in {Train, Test_original};
#   40 ASD / 31 CTRL; reference level = CTRL). Fits, in one R session and in this order:
#     1. sexadj : full ~ Batch + Sex + Diagnosis  vs  null ~ Batch + Sex
#     2. unadj  : full ~ Batch + Diagnosis        vs  null ~ Batch        (primary model)
#   Optional (--compadj, separate session in the original):
#     3. compadj: full ~ Batch + Sex + NLRproxy + Diagnosis vs null ~ Batch + Sex + NLRproxy
#   Settings: DRIMSeq 1.34.0; dmFilter(min_samps_gene_expr = 7, min_gene_expr = 10,
#   min_samps_feature_expr = 3, min_feature_expr = 10) -> 3,631 genes / 9,733 sites;
#   default dmPrecision (prec_subset = 0.1, trended moderation); gene- and feature-level
#   LRT; BH adjustment (DRIMSeq results()). Primary model: 81 genes adj_pvalue < 0.05
#   (2,217 genes with a finite P).
#
# Supports
#   Fig. 4B (blood DTU volcano; 81-gene signature) and fig. S4C (gene groups). The
#   sex-adjusted and composition-adjusted fits are sensitivity checks not shown in the paper. Read by
#   figures/fig4/fig5c_convergence.ipynb, analysis/blood/07-260615-within_family_systematic_apa.ipynb,
#   figures/fig4/fig5b_dtu_volcano.ipynb. Delta-psi / delta-length per gene
#   (gene_dpsi_curated.csv, gene_dlen_curated.csv) come from the companion script
#   02-261007-DTU_curated_genelevel_dpsi_dlen.py.
#
# Reproducibility
#   dmPrecision estimates the common precision on a random 10% subset of genes
#   (base::sample in the main process; BiocParallel 1.40 does not touch the global RNG).
#   The archived outputs are reproduced exactly only with set.seed(123) once at the top
#   and the sex-adjusted model fitted before the primary model in the same session
#   (i.e. the primary model uses the 2nd random draw). Fitting the primary model first
#   (1st draw) gives 80 instead of 81 hits. Do not change the seed or the model order.
#   compadj was written by a separate run (7 min later) whose RNG state is unknown:
#   --compadj runs it in a fresh session with set.seed(123) and is not bit-exact (see
#   Verification). Number of workers does not affect the results.
#
# Verification (2026-10-07, run with --outdir in a scratch dir)
#   dtu_unadj_gene/feature.csv and dtu_sexadj_gene/feature.csv: byte-identical to the
#   archived files (81 / 4 gene hits; 2,217 / 3,215 finite gene P; max |d log10 P| = 0).
#   dtu_compadj_*.csv: near only - same single gene hit, but per-gene P differ (max
#   |d log10 P| ~0.26; finite gene P 3,124 vs 3,112). Tried set.seed(123) with 0-3 prior
#   dmPrecision draws; none exact. NLRproxy itself was not reconstructed.
#
# Inputs
#   $M26_DATA_ROOT/3utr/Blood/all_till_20240724/utrcalling/
#       UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad
#       (136 samples x 26,142 sites / 8,338 genes; var utr_name = Ensembl gene id)
#   --compadj only: 2025_analyses/DTU_curated_sexadj/nlr_proxy_curated.csv
#       (per-sample UTR-derived neutrophil:lymphocyte proxy; its generating code was not
#       archived and is not reconstructed here - used as a fixed input)
#
# Outputs (default outdir = $M26_DATA_ROOT/3utr/Blood/2025_analyses/DTU_curated_sexadj)
#   dtu_sexadj_gene.csv, dtu_sexadj_feature.csv
#   dtu_unadj_gene.csv,  dtu_unadj_feature.csv
#   dtu_compadj_gene.csv, dtu_compadj_feature.csv   (only with --compadj)
#   Each table: DRIMSeq results() sorted by adj_pvalue (stable; NA last), write.csv(row.names = FALSE).
#   Existing files are never overwritten unless --overwrite is given.
#
# How to run (conda env apa-r; ~10 min for both models with 8 workers)
#   conda run -n apa-r Rscript 01-261007-DTU_curated_genelevel_unadj.R --outdir /path/to/scratch
#   conda run -n apa-r Rscript 01-261007-DTU_curated_genelevel_unadj.R --outdir /path --workers 8
#   conda run -n apa-r Rscript 01-261007-DTU_curated_genelevel_unadj.R --outdir /path --compadj
#   Options: --outdir DIR, --workers N (default 8), --compadj, --overwrite
#
# Provenance
#   Reconstructed on 2026-10-07 to reproduce DTU_curated_sexadj/dtu_unadj_gene.csv,
#   dtu_unadj_feature.csv, dtu_sexadj_gene.csv, dtu_sexadj_feature.csv (written together
#   2026-06-14 08:13) and, approximately, dtu_compadj_*.csv (2026-06-14 08:20); original
#   generating code was not archived. Data loading / DRIMSeq calls follow the
#   project's LOBO-CV script of the same period
#   (260614-LOBOCV_genelevel_curated_ridge_NOsexadj.R, not included).
# =============================================================================

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressMessages({
  library(DRIMSeq)
  library(dplyr)
})

# ---- arguments ---------------------------------------------------------------
parse_args <- function(args) {
  opt <- list(
    outdir = paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/DTU_curated_sexadj"),
    workers = 8L,
    compadj = FALSE,
    overwrite = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--outdir") {
      opt$outdir <- args[i + 1]
      i <- i + 1
    } else if (a == "--workers") {
      opt$workers <- as.integer(args[i + 1])
      i <- i + 1
    } else if (a == "--compadj") {
      opt$compadj <- TRUE
    } else if (a == "--overwrite") {
      opt$overwrite <- TRUE
    } else {
      stop("unknown argument: ", a)
    }
    i <- i + 1
  }
  opt
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

H5 <- paste0(
  paste0(DATA_ROOT, "/3utr/Blood/all_till_20240724/utrcalling/"),
  "UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad"
)
NLR_CSV <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/DTU_curated_sexadj/nlr_proxy_curated.csv")

set.seed(123) # must stay here, before any dmPrecision call (see header)

# ---- data: curated 71 samples, gene-level site counts ------------------------
adata <- anndata::read_h5ad(H5)
X <- adata$X
if (inherits(X, "sparseMatrix")) X <- as.matrix(X)
obs <- adata$obs
keep <- obs$status_split %in% c("Train", "Test_original")

meta <- data.frame(
  sample_id = gsub("-", ".", obs$sample_id[keep]),
  Diagnosis = relevel(factor(obs$Diagnosis[keep]), ref = "CTRL"),
  Batch = droplevels(factor(obs$Batch[keep])),
  Sex = factor(obs$Sex[keep])
)
cat(
  "curated n =", nrow(meta), "| ASD/CTRL:", sum(meta$Diagnosis == "ASD"), "/",
  sum(meta$Diagnosis == "CTRL"), "\n"
)

counts <- t(X[keep, , drop = FALSE])
colnames(counts) <- meta$sample_id
cdf <- data.frame(
  gene_id = adata$var$utr_name,
  feature_id = adata$var_names,
  counts,
  check.names = FALSE
)
rownames(cdf) <- NULL

d0 <- dmDSdata(counts = cdf, samples = meta)
d0 <- dmFilter(
  d0,
  min_samps_gene_expr = 7,
  min_gene_expr = 10,
  min_samps_feature_expr = 3,
  min_feature_expr = 10
)
cat("genes after dmFilter:", length(d0), "\n")

# ---- DRIMSeq fit + LRT -------------------------------------------------------
run_dtu <- function(d, full, null, smp, workers) {
  bp <- BiocParallel::MulticoreParam(workers = workers)
  design_full <- model.matrix(full, data = smp)
  design_null <- model.matrix(null, data = smp)
  d <- dmPrecision(d, design = design_full, BPPARAM = bp)
  d <- dmFit(d, design = design_full, BPPARAM = bp)
  d <- dmTest(d, design = design_null, BPPARAM = bp)
  list(
    gene = results(d, level = "gene") |> arrange(adj_pvalue),
    feature = results(d, level = "feature") |> arrange(adj_pvalue)
  )
}

write_res <- function(res, tag, outdir, overwrite) {
  for (lvl in c("gene", "feature")) {
    f <- file.path(outdir, paste0("dtu_", tag, "_", lvl, ".csv"))
    if (file.exists(f) && !overwrite) {
      cat("exists, not overwritten:", f, "\n")
      next
    }
    write.csv(res[[lvl]], f, row.names = FALSE)
    cat("wrote", f, "\n")
  }
  g <- res$gene
  cat(
    tag, "| genes tested:", nrow(g), "| finite P:", sum(!is.na(g$pvalue)),
    "| adj < 0.05:", sum(g$adj_pvalue < 0.05, na.rm = TRUE),
    "| adj < 0.1:", sum(g$adj_pvalue < 0.1, na.rm = TRUE), "\n"
  )
}

smp <- samples(d0)

if (!opt$compadj) {
  # order matters for exact reproduction: sexadj (1st RNG draw), then unadj (2nd draw)
  res_sex <- run_dtu(d0, ~ Batch + Sex + Diagnosis, ~ Batch + Sex, smp, opt$workers)
  write_res(res_sex, "sexadj", opt$outdir, opt$overwrite)

  res_unadj <- run_dtu(d0, ~ Batch + Diagnosis, ~Batch, smp, opt$workers)
  write_res(res_unadj, "unadj", opt$outdir, opt$overwrite)
} else {
  nlr <- read.csv(NLR_CSV)
  smp$NLRproxy <- nlr$NLRproxy[match(smp$sample_id, nlr$sample_id)]
  stopifnot(!anyNA(smp$NLRproxy))
  res_comp <- run_dtu(
    d0, ~ Batch + Sex + NLRproxy + Diagnosis, ~ Batch + Sex + NLRproxy, smp, opt$workers
  )
  write_res(res_comp, "compadj", opt$outdir, opt$overwrite)
}

sessionInfo() |> capture.output() |> head(12) |> cat(sep = "\n")
