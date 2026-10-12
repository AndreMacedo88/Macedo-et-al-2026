#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# 02-261012-gene_models.R
#
# Purpose
#   Per-gene test of diagnosis in the post-mortem cortex (Fig. 3C to 3E). For every
#   gene and cell type with a value in more than 10 samples, the per-sample mean
#   3'UTR length is modelled as length ~ Sex + Region + Diagnose and compared with
#   length ~ Sex + Region by a likelihood-ratio test (anova(..., test = "LRT")).
#   Genes whose model cannot be fitted (for example, all samples of one sex) get NA.
#
# Inputs (<indir>)
#   utr_per_cell_type_matrix.rda   mat_list, lm_tbl (01-261012-utr_lengths_per_cell_type.R)
#
# Outputs (<outdir>)
#   asd_utr_data.rda   lm_list (one named vector of P values per cell type), mat_list, lm_tbl
#
# How to run (conda env apa-r; about 5 min with 16 workers)
#   Rscript 02-261012-gene_models.R --indir DIR --outdir DIR [--workers 16]
#   Defaults: --indir and --outdir $M26_RESULTS_ROOT/cortex (or <repo>/results/cortex)
#
# Verification (2026-10-12)
#   On the matrices of 01-261012-utr_lengths_per_cell_type.R, lm_list matches the stored
#   asd_utr_data.rda of the 2020 analysis for all 244,074 tests to within 3e-8, except
#   one OPC region measured in 12 samples, which now gets P = 1 instead of NA.
#
# Provenance
#   Rewritten from human_asd_mapping.Rmd (2020 preprint of this work).
# =============================================================================

suppressMessages({
  library(parallel)
})

results_root <- Sys.getenv("M26_RESULTS_ROOT", "results")
opt <- list(indir = file.path(results_root, "cortex"), outdir = NULL, workers = "16")
args <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- args[i + 1]
}
if (is.null(opt$outdir)) opt$outdir <- opt$indir
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

load(file.path(opt$indir, "utr_per_cell_type_matrix.rda")) # mat_list, lm_tbl

gene_p <- function(x) {
  ok <- !is.na(x)
  if (sum(ok) <= 10) {
    return(NA_real_)
  }
  d <- data.frame(len = x, Diagnose = lm_tbl$Diagnose, Region = lm_tbl$Region, Sex = lm_tbl$Sex)
  fit1 <- tryCatch(lm(len ~ Sex + Region + Diagnose, data = d), error = function(e) NULL)
  if (is.null(fit1)) {
    return(NA_real_)
  }
  fit0 <- lm(len ~ Sex + Region, data = d)
  anova(fit1, fit0, test = "LRT")$`Pr(>Chi)`[2]
}

lm_list <- lapply(mat_list, function(m) {
  p <- unlist(mclapply(seq_len(ncol(m)), function(j) gene_p(m[, j]), mc.cores = as.integer(opt$workers)))
  setNames(p, colnames(m))
})

save(lm_list, mat_list, lm_tbl, file = file.path(opt$outdir, "asd_utr_data.rda"))
message(
  "P values for ", sum(vapply(lm_list, function(p) sum(!is.na(p)), numeric(1))),
  " gene x cell-type tests written to ", opt$outdir
)
