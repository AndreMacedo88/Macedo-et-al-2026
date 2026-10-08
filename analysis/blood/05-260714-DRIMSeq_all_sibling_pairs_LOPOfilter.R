#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# DRIMSeq DTU on all 9 in-house sibling pairs (18 samples): supp_8 (Fig. 4D, 4E) method
# (genewise precision, gene-level results, paired ~family+Diagnosis) re-run at the
# LOPO expression filter (and, for context, across the full filter sweep).
# -----------------------------------------------------------------------------
# Background. An earlier version of supp_8 (figures/fig4/supp_8_within_family_drimseq.ipynb)
# reported 8 genes adj-p<0.05 from an unarchived script with a lenient expression filter.
# The filter sweep below shows that the
# supp_8 method combined with the tighter LOPO filter is a valid independent-filtering
# refinement that recovers more genes after BH (8 -> ~20). This script regenerates that
# result from scratch and stores it.
#
# Main output = the LOPO-filter run:
#   drimseq_allpairs_genelevel_LOPOfilter.csv   (gene-level, used by supp_8; ~20 adj<0.05)
#   drimseq_allpairs_featurelevel_LOPOfilter.csv
# Context output = the same method across 5 filters (very_lenient .. tightest), one
#   gene-level table per filter (genelevel_<filter>.csv) + a long significant-gene table
#   (sig_genes_by_filter.csv) that the overlap notebook
#   260714-filtersweep_vs_LOPO_overlap_analysis.ipynb (not included) reads.
#
# Method (matches supp_8; the only non-default is prec_moderation="none", which is
#   required because the default "trended" loess fails on this small 18-sample design;
#   this setting reproduces the 8 genes of the lenient-filter version, see the
#   very_lenient row of the filter sweep):
#     dmDSdata -> dmFilter -> dmPrecision(genewise, moderation="none") ->
#     dmFit(~family+Diagnosis) -> dmTest(null ~family) -> results(level="gene"/"feature").
#   Not common precision, not feature-then-gene: this is the stable, gene-level,
#   FDR-correct within-family DTU readout (distinct from the LOPO classifier's
#   common-precision feature-level selection).
#
# Env: conda apa-r. Usage: Rscript 05-260714-DRIMSeq_all_sibling_pairs_LOPOfilter.R
# =============================================================================
DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressMessages({ library(DRIMSeq) })
Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
BASE <- paste0(DATA_ROOT, "/3utr/Blood/all_till_20240724/utrcalling/")
COUNTS_H5 <- paste0(BASE, "UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad")
HERE <- paste0(DATA_ROOT, "/notebooks/blood_samples/Classifications/2025/LOPO_siblings")
SHEET <- file.path(HERE, "samplesheet_inhouse_9pairs.csv")
OUTDIR <- paste0(DATA_ROOT, "/3utr/Blood/2025_analyses/DRIMSeq_allpairs_genewise")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
SEED <- 123L; CORES <- 8L
BPPARAM <- BiocParallel::MulticoreParam(workers = CORES)

# Expression filters (min_samps_gene_expr / min_gene_expr / min_samps_feature_expr / min_feature_expr)
LOPO_FILT <- list(min_samps_gene_expr = 6, min_gene_expr = 10, min_samps_feature_expr = 3, min_feature_expr = 10)
FILTERS <- list(
  very_lenient = list(min_samps_gene_expr = 3, min_gene_expr = 10, min_samps_feature_expr = 2, min_feature_expr = 5),
  lenient      = list(min_samps_gene_expr = 3, min_gene_expr = 10, min_samps_feature_expr = 3, min_feature_expr = 10),
  LOPO_tight   = LOPO_FILT,
  tighter      = list(min_samps_gene_expr = 9, min_gene_expr = 10, min_samps_feature_expr = 6, min_feature_expr = 10),
  tightest     = list(min_samps_gene_expr = 12, min_gene_expr = 20, min_samps_feature_expr = 9, min_feature_expr = 20)
)

read_h5ad_dense <- function(path) {
  ad_py <- reticulate::import("anndata", convert = FALSE); sp <- reticulate::import("scipy.sparse", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE); A <- ad_py$read_h5ad(path); Xpy <- A$X
  if (isTRUE(reticulate::py_to_r(sp$issparse(Xpy)))) Xpy <- Xpy$toarray()
  X <- reticulate::py_to_r(np$asarray(Xpy, dtype = "float64"))
  rownames(X) <- as.character(reticulate::py_to_r(A$obs_names$tolist()))
  colnames(X) <- as.character(reticulate::py_to_r(A$var_names$tolist()))
  var_r <- reticulate::py_to_r(A$var)
  list(Xfs = t(X), gene = as.character(var_r$utr_name), sym = as.character(var_r$gene_name),
       features = colnames(X), samples = rownames(X))
}
rc <- read_h5ad_dense(COUNTS_H5)
gene_id <- rc$gene; feature_id <- rc$features
# ENSG group -> gene symbol map (utr_name == name; gene_name = symbol)
ensg2sym <- tapply(rc$sym, rc$gene, function(v) v[1])

sheet <- read.csv(SHEET, stringsAsFactors = FALSE)
ids <- sheet$sample_id
stopifnot(all(ids %in% rc$samples))
counts_mat <- round(rc$Xfs[, ids, drop = FALSE])
meta_all <- data.frame(sample_id = sheet$sample_id,
                       Diagnosis = factor(sheet$Diagnosis, levels = c("CTRL", "ASD")),
                       family = factor(sheet$family), stringsAsFactors = FALSE)
rownames(meta_all) <- meta_all$sample_id
message("all sibling pairs: ", nrow(meta_all), " samples | ASD=", sum(meta_all$Diagnosis == "ASD"),
        " CTRL=", sum(meta_all$Diagnosis == "CTRL"), " | families=", nlevels(meta_all$family))

# genewise, gene- and feature-level, at a given filter (supp_8 method)
fit_supp8 <- function(filt) {
  set.seed(SEED)  # reproducibility for dmPrecision's common-precision init (10% subset)
  df <- data.frame(gene_id = gene_id, feature_id = feature_id, counts_mat, check.names = FALSE)
  d <- dmDSdata(counts = df, samples = meta_all)
  d <- dmFilter(d, min_samps_gene_expr = filt$min_samps_gene_expr, min_gene_expr = filt$min_gene_expr,
                min_samps_feature_expr = filt$min_samps_feature_expr, min_feature_expr = filt$min_feature_expr)
  n_genes_filt <- length(unique(counts(d)$gene_id)); n_feats_filt <- nrow(counts(d))
  design_full <- model.matrix(~ family + Diagnosis, data = samples(d))
  design_null <- model.matrix(~ family, data = samples(d))
  d <- dmPrecision(d, design = design_full, prec_moderation = "none", BPPARAM = BPPARAM)   # genewise
  d <- dmFit(d, design = design_full, BPPARAM = BPPARAM)
  d <- dmTest(d, design = design_null, BPPARAM = BPPARAM)
  rG <- results(d, level = "gene");    rG <- rG[!is.na(rG$pvalue), ]
  rF <- results(d, level = "feature"); rF <- rF[!is.na(rF$pvalue), ]
  rG <- rG[order(rG$pvalue), ]; rF <- rF[order(rF$pvalue), ]
  rG$symbol <- unname(ensg2sym[as.character(rG$gene_id)])
  rF$symbol <- unname(ensg2sym[as.character(rF$gene_id)])
  list(gene = rG, feat = rF, n_genes_filt = n_genes_filt, n_feats_filt = n_feats_filt)
}

sig_long <- list(); summ <- list()
for (nm in names(FILTERS)) {
  t0 <- Sys.time()
  r <- fit_supp8(FILTERS[[nm]])
  # reorder gene-level columns: gene_id, symbol, lr, df, pvalue, adj_pvalue
  g <- r$gene[, c("gene_id", "symbol", "lr", "df", "pvalue", "adj_pvalue")]
  write.csv(g, file.path(OUTDIR, paste0("genelevel_", nm, ".csv")), row.names = FALSE)
  s05 <- g[!is.na(g$adj_pvalue) & g$adj_pvalue < 0.05, ]
  sig_long[[nm]] <- data.frame(filter = nm, gene_id = g$gene_id, symbol = g$symbol,
                               pvalue = g$pvalue, adj_pvalue = g$adj_pvalue,
                               sig05 = !is.na(g$adj_pvalue) & g$adj_pvalue < 0.05,
                               sig10 = !is.na(g$adj_pvalue) & g$adj_pvalue < 0.10,
                               stringsAsFactors = FALSE)
  summ[[nm]] <- data.frame(filter = nm, n_genes_filt = r$n_genes_filt, n_feats_filt = r$n_feats_filt,
                           gene_tested = nrow(g),
                           gene_sig05 = sum(g$adj_pvalue < 0.05, na.rm = TRUE),
                           gene_sig10 = sum(g$adj_pvalue < 0.10, na.rm = TRUE))
  # Main result: the LOPO-filter run (supp_8 method at the tighter filter)
  if (nm == "LOPO_tight") {
    write.csv(g, file.path(OUTDIR, "drimseq_allpairs_genelevel_LOPOfilter.csv"), row.names = FALSE)
    write.csv(r$feat[, c("gene_id", "symbol", "feature_id", "lr", "df", "pvalue", "adj_pvalue")],
              file.path(OUTDIR, "drimseq_allpairs_featurelevel_LOPOfilter.csv"), row.names = FALSE)
  }
  message(sprintf("  [%-12s] genes_tested=%4d  gene adj<.05=%2d  adj<.10=%2d  (%.0fs)",
                  nm, nrow(g), sum(g$adj_pvalue < 0.05, na.rm = TRUE),
                  sum(g$adj_pvalue < 0.10, na.rm = TRUE), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
write.csv(do.call(rbind, sig_long), file.path(OUTDIR, "sig_genes_by_filter.csv"), row.names = FALSE)
write.csv(do.call(rbind, summ), file.path(OUTDIR, "filter_sweep_summary.csv"), row.names = FALSE)

hl <- summ[["LOPO_tight"]]
writeLines(c(
  "DRIMSeq DTU on all 9 in-house sibling pairs (18 samples) - supp_8 method (genewise, gene-level)",
  paste("counts     :", COUNTS_H5),
  paste("samplesheet:", SHEET),
  "design     : full ~ family + Diagnosis ; null ~ family",
  "precision  : genewise (prec_moderation = none)",
  "HEADLINE FILTER (LOPO): min_samps_gene_expr=6, min_gene_expr=10, min_samps_feature_expr=3, min_feature_expr=10",
  paste("HEADLINE genes tested :", hl$gene_tested),
  paste("HEADLINE gene adj<0.05:", hl$gene_sig05, " | adj<0.10:", hl$gene_sig10),
  "(supp_8 original at its lenient filter = 8 adj<0.05 / 20 adj<0.10; reproduced by the very_lenient row)",
  paste("date       :", as.character(Sys.time()))
), file.path(OUTDIR, "run_config.txt"))

message("HEADLINE (LOPO filter, supp_8 re-run): ", hl$gene_sig05, " genes adj<0.05, ", hl$gene_sig10, " adj<0.10")
message("Outputs in: ", OUTDIR)
