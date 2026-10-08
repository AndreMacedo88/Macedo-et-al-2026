#!/usr/bin/env Rscript

# =============================================================================
# Pseudotime ORA on the chi-square + permutation intersection foreground
# (Reactome and DisGeNET), Knoblich human brain organoids, both lineages.
#
# Supports: Fig. 2E (Reactome / DisGeNET ORA of developmentally regulated
#   3'UTR genes, INH; EXC counterpart for supplement / text).
#
# Rewritten on 2026-10-07; reproduces the stored
#   ORA/chisq_perm_intersection/ora_{reactome,dgn}_{inhibitory,excitatory}_neglog10pval_chisqperm.rds.
#   The gene filtering, Entrez
#   mapping and enrichment calls follow 14_260506-knoblich_enrichment_analyses.ipynb
#   (chi-square-only ORA); only the foreground definition differs.
#
# Foreground / background (per lineage):
#   1. composite_ranking_optimized2_<lin>.csv, filtered exactly as in
#      14_260506-knoblich_enrichment_analyses.ipynb: converged genes, no
#      comma-separated gene_id, leading '|' stripped from gene_name, rows whose
#      gene_name still contains '|' removed.
#   2. padj_chisq is recomputed (BH) on p_value_chisq over this filtered set
#      (not the padj / padj_chisq columns of the composite CSV; using those
#      makes the INH foreground 12 genes too large).
#   3. The within-sample permutation p-value (p_value_perm, joined by gene_id)
#      is BH-adjusted over the same filtered set (not the padj column of the
#      permutation CSV; using that makes the EXC foreground too large).
#   4. Foreground = Entrez IDs of genes with padj_chisq < 0.05 and
#      padj_perm < 0.05 (INH 169, EXC 565 Entrez IDs).
#      Background = all Entrez IDs of the filtered set (INH 8,341; EXC 9,655).
#   Entrez mapping: org.Hs.eg.db, ENSEMBL -> ENTREZID, multiVals = "first".
#   Enrichment: ReactomePA::enrichPathway and DOSE::enrichDGN, BH,
#   minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1 (all terms kept).
#
# Inputs (KB = $M26_DATA_ROOT/3utr/organoids/knoblich/utrcalling):
#   KB/composite_ranking/composite_ranking_optimized2_<lin>.csv
#   KB/dirmult_results/<lin>/permutation_batch_number/
#     dirmult_regression_permutation_optimized2_joint_df5_pt-dpt_pseudotime_<lin>_cov-batch_number.csv
#
# Outputs (default KB/ORA/chisq_perm_intersection/):
#   ora_reactome_<lin>_neglog10pval_chisqperm.rds  (enrichResult)
#   ora_dgn_<lin>_neglog10pval_chisqperm.rds       (enrichResult)
#
# Run (conda env apa-r):
#   conda run -n apa-r Rscript 15_261007-pseudotime_ora_chisq_perm_intersection.R \
#     [--outdir <dir>] [--lineages inhibitory,excitatory]
#
# Verified 2026-10-07 (clusterProfiler 4.14.0, ReactomePA 1.50.0, DOSE 4.0.0,
# org.Hs.eg.db 3.20.0, reactome.db 1.89.0): @gene, @universe and @result
# (ID, pvalue, p.adjust, Count, geneID) identical to the archived RDS files.
# =============================================================================

DATA_ROOT <- Sys.getenv("M26_DATA_ROOT")  # data directory (see README)
if (!nzchar(DATA_ROOT)) stop("Set M26_DATA_ROOT (see README)")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(DOSE)
  library(dplyr)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(readr)
  library(stringr)
})

KB <- paste0(DATA_ROOT, "/3utr/organoids/knoblich/utrcalling")
PADJ_CUT <- 0.05
METRIC <- "neglog10pval_chisqperm"

parse_args <- function(args) {
  opts <- list(
    outdir = file.path(KB, "ORA", "chisq_perm_intersection"),
    lineages = c("inhibitory", "excitatory")
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    val <- args[[i + 1]]
    if (key == "--outdir") {
      opts$outdir <- val
    } else if (key == "--lineages") {
      opts$lineages <- strsplit(val, ",")[[1]]
    } else {
      stop(sprintf("Unknown argument: %s", key))
    }
    i <- i + 2
  }
  opts
}

# Filtering identical to load_and_filter() in 14_260506-knoblich_enrichment_analyses.ipynb
load_and_filter <- function(lin) {
  path <- file.path(KB, "composite_ranking", sprintf("composite_ranking_optimized2_%s.csv", lin))
  dm <- read_csv(path, show_col_types = FALSE) |>
    filter(converged == TRUE) |>
    filter(!grepl(",", gene_id)) |>
    mutate(gene_name = str_trim(sub("^\\|+", "", gene_name))) |>
    filter(!grepl("\\|", gene_name) | is.na(gene_name)) |>
    mutate(padj_chisq = p.adjust(p_value_chisq, method = "BH"))
  dm$entrezid <- mapIds(
    org.Hs.eg.db,
    keys = dm$gene_id, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"
  )
  dm
}

add_permutation <- function(dm, lin) {
  path <- file.path(
    KB, "dirmult_results", lin, "permutation_batch_number",
    sprintf("dirmult_regression_permutation_optimized2_joint_df5_pt-dpt_pseudotime_%s_cov-batch_number.csv", lin)
  )
  pm <- read_csv(path, show_col_types = FALSE)
  dm |>
    mutate(
      p_value_perm = pm$p_value_perm[match(gene_id, pm$utr_name)],
      padj_perm = p.adjust(p_value_perm, method = "BH")
    )
}

run_lineage <- function(lin, outdir) {
  dm <- load_and_filter(lin) |> add_permutation(lin)

  fg <- dm |>
    filter(padj_chisq < PADJ_CUT, padj_perm < PADJ_CUT) |>
    pull(entrezid) |>
    na.omit() |>
    unique() |>
    as.character()
  bg <- dm$entrezid[!is.na(dm$entrezid)] |>
    unique() |>
    as.character()
  cat(sprintf(
    "[%s] filtered genes %d | chisq&perm sig genes %d | foreground %d Entrez | background %d Entrez\n",
    lin, nrow(dm), sum(dm$padj_chisq < PADJ_CUT & dm$padj_perm < PADJ_CUT, na.rm = TRUE),
    length(fg), length(bg)
  ))

  ora_reactome <- enrichPathway(
    gene = fg, universe = bg, organism = "human",
    pvalueCutoff = 1, pAdjustMethod = "BH",
    minGSSize = 15, maxGSSize = 500
  )
  ora_dgn <- enrichDGN(
    gene = fg, universe = bg,
    pvalueCutoff = 1, pAdjustMethod = "BH",
    minGSSize = 15, maxGSSize = 500
  )

  for (db in c("reactome", "dgn")) {
    obj <- if (db == "reactome") ora_reactome else ora_dgn
    out <- file.path(outdir, sprintf("ora_%s_%s_%s.rds", db, lin, METRIC))
    saveRDS(obj, out)
    n_annot <- sum(obj@gene %in% unique(unlist(obj@geneSets)))
    cat(sprintf(
      "[%s] %s: %d terms tested, %d with p.adjust < 0.05, %d / %d foreground genes annotated -> %s\n",
      lin, db, nrow(obj@result), sum(obj@result$p.adjust < 0.05), n_annot, length(obj@gene), out
    ))
  }
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)
for (lin in opts$lineages) run_lineage(lin, opts$outdir)
sessionInfo() |> capture.output() |> head(3) |> cat(sep = "\n")
