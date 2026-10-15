# Table S1: Gene Ontology over-representation behind Fig. 1F
#
# Exports every GO biological-process term tested in the over-representation analysis of
# Fig. 1F (foreground: genes significant under both the chi-square LRT and the permutation
# test; background: tested genes with a GO annotation), checks it against the values given in
# the Results, and writes one CSV row per term.
#
# Input:
#   DirMult_optimized/ORA/chisq_perm_intersection/ora_go_neglog10pval_chisqperm.rds
#     (written by analysis/mouse_svz/09_261007-optimized2_ora_chisq_perm_intersection.R)
# Output (default: $M26_RESULTS_ROOT/tables, or <repo>/results/tables):
#   table_s1_go_enrichment.csv (converted to the Excel table by table_s1_xlsx.py)
#
# How to run:
#   conda run -n apa-r Rscript table_s1_go_enrichment.R [--outdir DIR]

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT") # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
})

ora_path <- file.path(
  DATA_ROOT, "3utr/wildtypes_analysis/2026_analyses/utrcalling",
  "DirMult_optimized/ORA/chisq_perm_intersection/ora_go_neglog10pval_chisqperm.rds"
)
args <- commandArgs(trailingOnly = TRUE)
outdir <- file.path(Sys.getenv("M26_RESULTS_ROOT", "results"), "tables")
if ("--outdir" %in% args) outdir <- args[which(args == "--outdir") + 1]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(outdir, "table_s1_go_enrichment.csv")

# terms drawn in Fig. 1F (figures/fig1/fig1f_bp_enrichment)
fig1f_terms <- c(
  "regulation of RNA splicing",
  "mRNA processing",
  "cognition",
  "learning or memory",
  "cytoplasmic translation",
  "apoptotic mitochondrial changes"
)

ora <- readRDS(ora_path)
res <- ora@result
stopifnot(ora@ontology == "BP", length(ora@gene) == 1746, length(ora@universe) == 7214)

split_ratio <- function(x, part) as.integer(sapply(strsplit(x, "/"), `[`, part))

tab <- res |>
  mutate(
    fg_in_term = split_ratio(GeneRatio, 1),
    fg_annotated = split_ratio(GeneRatio, 2),
    bg_in_term = split_ratio(BgRatio, 1),
    bg_annotated = split_ratio(BgRatio, 2),
    fold_enrichment = (fg_in_term / fg_annotated) / (bg_in_term / bg_annotated),
    shown_in_fig1F = ifelse(Description %in% fig1f_terms, "yes", ""),
    genes = sapply(strsplit(geneID, "/"), function(g) paste(sort(g), collapse = ", "))
  ) |>
  arrange(pvalue, Description) |>
  transmute(
    go_id = ID,
    term = Description,
    shown_in_fig1F,
    fg_in_term,
    fg_annotated,
    bg_in_term,
    bg_annotated,
    fold_enrichment,
    p_value = pvalue,
    p_adjusted_BH = p.adjust,
    genes
  )

# consistency with the stored object and with the values given in the Results
stopifnot(
  nrow(tab) == 1436,
  all(tab$fg_in_term == res$Count[match(tab$go_id, res$ID)]),
  all(fig1f_terms %in% tab$term),
  all(tab$fg_annotated == 1699),
  all(tab$bg_annotated == 7214)
)
quoted <- tibble::tribble(
  ~term, ~padj, ~n,
  "regulation of RNA splicing", 1.1e-4, 53,
  "mRNA processing", 8.2e-4, 106,
  "cognition", 0.0124, 66,
  "learning or memory", 0.040, 58,
  "cytoplasmic translation", 4.2e-3, 28
)
chk <- quoted |>
  left_join(tab, by = "term") |>
  mutate(ok = signif(p_adjusted_BH, 2) == signif(padj, 2) & fg_in_term == n)
print(chk |> select(term, padj, p_adjusted_BH, n, fg_in_term, ok))
stopifnot(all(chk$ok))

write.csv(tab, out_path, row.names = FALSE)
cat(
  "terms:", nrow(tab), "| adjusted P < 0.05:", sum(tab$p_adjusted_BH < 0.05),
  "| shown in Fig. 1F:", sum(tab$shown_in_fig1F == "yes"), "\nwritten:", out_path, "\n"
)
