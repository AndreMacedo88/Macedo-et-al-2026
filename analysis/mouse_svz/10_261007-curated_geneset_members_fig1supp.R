# Curated functional gene-set membership for fig. S3A (mouse adult SVZ, optimized method 2)
#
# Purpose:
#   Define the three curated functional gene sets whose mean 3'UTR-length index W(t) is
#   plotted along pseudotime in fig. S3A (figures/fig1/supp_4_exemplars_genesets.ipynb,
#   cell "S4b"). Each set is the union of the ORA foreground genes (geneID column) of the
#   curated GO BP terms shown in the Fig. 1 biological-process panel
#   (Fig. 1F, figures/fig1/fig1f_bp_enrichment.ipynb, `curated_bio`), taken from the base chi-square-only GO
#   ORA (foreground = padj_chisq < 0.05). Genes keep their order of first appearance
#   (terms in the order listed below; within a term, the ORA geneID order = rank order).
#     RNA processing   (119 genes): regulation of RNA splicing; RNA splicing;
#                                   regulation of mRNA splicing, via spliceosome;
#                                   alternative mRNA splicing, via spliceosome; mRNA processing
#     Cognition/neural (154 genes): cognition; learning or memory; memory; nervous system process
#     Proteostasis     (156 genes): ubiquitin-dependent protein catabolic process;
#                                   cytoplasmic translation
#
# Supports: fig. S3A (functional gene-set W(t) trajectories).
#
# Input:
#   DirMult_optimized/ORA/ora_go_neglog10pval_chisq_optimized2.rds
#     (written by 08_260411-optimized2_enrichment_analyses.ipynb)
#
# Output (default: DirMult_optimized/ORA/):
#   curated_geneset_members_fig1supp.csv (columns: category, gene_name)
#
# How to run:
#   conda run -n apa-r Rscript 10_261007-curated_geneset_members_fig1supp.R [--outdir DIR] [--overwrite]
#
# Reconstructed on 2026-10-07 to reproduce DirMult_optimized/ORA/curated_geneset_members_fig1supp.csv;
# original generating code was not archived. The rule above reproduces the archived file
# exactly (same members, same order).

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT") # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
})

ora_dir <- paste0(DATA_ROOT, "/3utr/wildtypes_analysis/2026_analyses/utrcalling/DirMult_optimized/ORA")
ora_go_path <- file.path(ora_dir, "ora_go_neglog10pval_chisq_optimized2.rds")

args <- commandArgs(trailingOnly = TRUE)
outdir <- ora_dir
if ("--outdir" %in% args) outdir <- args[which(args == "--outdir") + 1]
out_file <- file.path(outdir, "curated_geneset_members_fig1supp.csv")
if (file.exists(out_file) && !("--overwrite" %in% args)) {
  stop(out_file, " exists; pass --overwrite to replace it")
}

curated_terms <- list(
  "RNA processing" = c(
    "regulation of RNA splicing",
    "RNA splicing",
    "regulation of mRNA splicing, via spliceosome",
    "alternative mRNA splicing, via spliceosome",
    "mRNA processing"
  ),
  "Cognition/neural" = c(
    "cognition",
    "learning or memory",
    "memory",
    "nervous system process"
  ),
  "Proteostasis" = c(
    "ubiquitin-dependent protein catabolic process",
    "cytoplasmic translation"
  )
)

ora_go <- readRDS(ora_go_path)@result
missing_terms <- setdiff(unlist(curated_terms), ora_go$Description)
if (length(missing_terms) > 0) stop("terms not found in GO ORA: ", paste(missing_terms, collapse = "; "))

members <- curated_terms |>
  imap(\(terms, category) {
    genes <- ora_go$geneID[match(terms, ora_go$Description)] |>
      strsplit("/") |>
      unlist() |>
      unique()
    tibble(category = category, gene_name = genes)
  }) |>
  list_rbind()

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(members, out_file, row.names = FALSE)
print(count(members, category))
cat("Wrote", out_file, "\n")
