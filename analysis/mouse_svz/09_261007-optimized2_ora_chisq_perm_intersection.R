# ORA on the chi-square / permutation intersection foreground (mouse adult SVZ, optimized method 2)
#
# Purpose:
#   Over-representation analysis (GO BP, Reactome, KEGG, MSigDB Hallmark) of genes whose
#   3'UTR isoform usage changes along pseudotime, using as foreground the genes that are
#   significant at BH < 0.05 under both the Dirichlet-multinomial chi-square LRT and the
#   within-sample permutation test. Same settings as the base pipeline
#   (08_260411-optimized2_enrichment_analyses.ipynb, chi-square-only ORA), except that the
#   foreground is the intersection and all terms are kept (pvalueCutoff = 1).
#
# Supports: Fig. 1F (GO BP) and fig. S1J (Reactome, KEGG, Hallmark).
#
# Inputs:
#   DirMult_optimized/dirmult_regression_optimized2_joint_df5.csv
#     (chi-square LRT; BH recomputed over converged, single-gene entries; also the universe)
#   DirMult_optimized/permutation_batch_number/
#     dirmult_regression_permutation_optimized2_joint_df5_pt-dpt_pseudotime_cov-batch_number.csv
#     (permutation test; column padj as written by 03_260526-permutation_optimized_2.R)
#   $M26_DATA_ROOT/common_files/annotations/mouse_to_human_genes.txt (Ensembl orthologs)
#   org.Mm.eg.db, org.Hs.eg.db, ReactomePA (reactome.db), msigdbr, KEGG rest (online)
#
# Outputs (default: DirMult_optimized/ORA/chisq_perm_intersection/):
#   ora_go_neglog10pval_chisqperm.rds
#   ora_reactome_neglog10pval_chisqperm.rds
#   ora_kegg_neglog10pval_chisqperm.rds
#   ora_hallmark_neglog10pval_chisqperm.rds
#
# How to run:
#   conda run -n apa-r Rscript 09_261007-optimized2_ora_chisq_perm_intersection.R [--outdir DIR] [--overwrite]
#
# Notes:
#   KEGG is queried online (enrichKEGG), so results depend on the KEGG release at run time
#   (2026-10-07 rerun: GO, Reactome, Hallmark identical to the archived rds; KEGG universe
#   3261 vs 3243, one extra term, max |delta p.adjust| 0.04, same 6 terms at FDR < 0.05).
#
# Reconstructed on 2026-10-07 to reproduce
# DirMult_optimized/ORA/chisq_perm_intersection/ora_{go,reactome,kegg,hallmark}_neglog10pval_chisqperm.rds;
# original generating code was not archived.

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ReactomePA)
  library(msigdbr)
})

# ---- paths and arguments ----
dmo_dir <- paste0(DATA_ROOT, "/3utr/wildtypes_analysis/2026_analyses/utrcalling/DirMult_optimized")
chisq_path <- file.path(dmo_dir, "dirmult_regression_optimized2_joint_df5.csv")
perm_path <- file.path(
  dmo_dir, "permutation_batch_number",
  "dirmult_regression_permutation_optimized2_joint_df5_pt-dpt_pseudotime_cov-batch_number.csv"
)
ortholog_path <- paste0(DATA_ROOT, "/common_files/annotations/mouse_to_human_genes.txt")

args <- commandArgs(trailingOnly = TRUE)
outdir <- file.path(dmo_dir, "ORA", "chisq_perm_intersection")
if ("--outdir" %in% args) outdir <- args[which(args == "--outdir") + 1]
metric <- "neglog10pval_chisqperm"
alpha <- 0.05
existing <- file.exists(file.path(outdir, paste0("ora_", c("go", "reactome", "kegg", "hallmark"), "_", metric, ".rds")))
if (any(existing) && !("--overwrite" %in% args)) {
  stop("outputs already exist in ", outdir, "; pass --overwrite to replace them")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- chi-square LRT table (identical preprocessing to the base notebook) ----
dm2 <- readr::read_csv(chisq_path, show_col_types = FALSE) |>
  filter(converged == TRUE) |>
  filter(!grepl(",", utr_name))
dm2$padj_chisq <- p.adjust(dm2$p_value_chisq, method = "BH")
min_p <- min(dm2$p_value_chisq[dm2$p_value_chisq > 0])
dm2$p_value_chisq[dm2$p_value_chisq == 0] <- min_p
dm2 <- dm2 |>
  mutate(rank_score = -log10(p_value_chisq)) |>
  arrange(desc(rank_score))

# ---- permutation test: join its BH-adjusted p-value ----
perm <- readr::read_csv(perm_path, show_col_types = FALSE) |>
  dplyr::select(utr_name, padj_perm = padj)
dm2 <- dm2 |> left_join(perm, by = "utr_name")
dm2$sig_both <- dm2$padj_chisq < alpha & !is.na(dm2$padj_perm) & dm2$padj_perm < alpha
cat("Converged genes:", nrow(dm2), "| chi-sq sig:", sum(dm2$padj_chisq < alpha),
    "| chi-sq AND perm sig:", sum(dm2$sig_both), "\n")

# ---- gene ID mappings (identical to the base notebook) ----
mouse_to_human_genes <- readr::read_csv(ortholog_path, show_col_types = FALSE)
dm2$gene_symbol <- mapIds(org.Mm.eg.db, dm2$utr_name, "SYMBOL", "ENSEMBL")
one2one <- mouse_to_human_genes |>
  filter(grepl("ortholog_one2one|ortholog_one2many", `Human homology type`))
one2one_in_data <- one2one[one2one$`Gene stable ID` %in% dm2$utr_name, ] |>
  arrange(
    desc(`Human orthology confidence [0 low, 1 high]`),
    desc(`Human Whole-genome alignment coverage`)
  ) |>
  filter(!duplicated(`Gene stable ID`))
human_ids <- one2one_in_data[match(dm2$utr_name, one2one_in_data$`Gene stable ID`), ]
dm2$entrezid <- as.numeric(as.character(
  mapIds(org.Hs.eg.db, human_ids$`Human gene stable ID`, "ENTREZID", "ENSEMBL")
))
dm2$mouse_entrezid <- mapIds(org.Mm.eg.db, dm2$utr_name, "ENTREZID", "ENSEMBL")

# ---- foreground / universe vectors ----
# Foreground genes are taken in Ensembl ID order (the order of both input CSVs), the universe
# in rank order (as in the base notebook). Order only affects the order of @gene and of the
# geneID strings, not the statistics; it is kept so the archived rds files are reproduced exactly.
fg_bg <- function(col) {
  fg <- dm2 |>
    filter(sig_both) |>
    arrange(utr_name) |>
    pull({{ col }}) |>
    na.omit() |>
    unique() |>
    as.character()
  bg <- dm2 |> pull({{ col }}) |> na.omit() |> unique() |> as.character()
  list(fg = fg, bg = bg)
}
sym <- fg_bg(gene_symbol)
hum <- fg_bg(entrezid)
kegg <- fg_bg(mouse_entrezid)
cat("Foreground symbols:", length(sym$fg), "of", length(sym$bg), "\n")

# ---- ORA (pvalueCutoff = 1 keeps every tested term) ----
ora_go <- enrichGO(
  gene = sym$fg, universe = sym$bg,
  OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
  pvalueCutoff = 1, qvalueCutoff = 0.2, pAdjustMethod = "BH",
  minGSSize = 50, maxGSSize = 500
)

ora_reactome <- enrichPathway(
  gene = hum$fg, universe = hum$bg, organism = "human",
  pvalueCutoff = 1, qvalueCutoff = 0.2, pAdjustMethod = "BH",
  minGSSize = 20, maxGSSize = 500
)

ora_kegg <- enrichKEGG(
  gene = kegg$fg, universe = kegg$bg, organism = "mmu",
  pvalueCutoff = 1, qvalueCutoff = 0.2, pAdjustMethod = "BH",
  minGSSize = 50, maxGSSize = 500
)
ora_kegg@result$Description <- gsub(
  " - Mus musculus \\(house mouse\\)$", "", ora_kegg@result$Description
)

hallmark_t2g <- msigdbr(species = "mouse", collection = "H") |>
  dplyr::select(gs_name, gene_symbol) |>
  dplyr::distinct()
hallmark_t2n <- hallmark_t2g |>
  dplyr::select(gs_name) |>
  dplyr::distinct() |>
  dplyr::mutate(
    gs_description = gs_name |>
      str_remove("^HALLMARK_") |>
      str_replace_all("_", " ") |>
      str_to_title()
  )
ora_hallmark <- enricher(
  gene = sym$fg, universe = sym$bg,
  TERM2GENE = hallmark_t2g, TERM2NAME = hallmark_t2n,
  pvalueCutoff = 1, qvalueCutoff = 0.2, pAdjustMethod = "BH",
  minGSSize = 20, maxGSSize = 500
)

# ---- save ----
results <- list(go = ora_go, reactome = ora_reactome, kegg = ora_kegg, hallmark = ora_hallmark)
for (db in names(results)) {
  res <- results[[db]]
  out_file <- file.path(outdir, paste0("ora_", db, "_", metric, ".rds"))
  saveRDS(res, out_file)
  cat(sprintf(
    "%-9s terms: %4d | FDR<0.05: %3d | -> %s\n",
    db, nrow(res@result), sum(res@result$p.adjust < alpha), out_file
  ))
}
