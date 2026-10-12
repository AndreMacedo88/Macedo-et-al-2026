#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# 03-261012-gsea_l23.R
#
# Purpose
#   Reactome gene-set enrichment for layer 2/3 excitatory neurons (Fig. 3E), as run for
#   the 2020 preprint of this work. Every 3'UTR region is ranked by the log ratio of its
#   mean length in ASD and control samples (one region per gene for most genes; genes with
#   several separate 3'UTR regions appear several times). A gene enters the pathways through
#   its highest-ranked region and its other regions stay in the ranking as non-members, which
#   is how the fgsea version of 2020 treated repeated gene identifiers. Simple permutation
#   test (fgseaSimple, 10,000 permutations, gene sets of 10 to 500 genes), Benjamini-Hochberg
#   adjustment; pathways with an adjusted P below 0.1 are kept.
#
# Inputs (<indir>)
#   asd_utr_data.rda   mat_list, lm_tbl (02-261012-gene_models.R)
#   reactome.db        Reactome pathways (Bioconductor annotation package)
#
# Outputs (<outdir>)
#   asd_pathway.rda           path_l23_tbl (pathways with adjusted P < 0.1)
#   gsea_l23_reactome.tsv     all tested pathways
#
# How to run (conda env apa-r; about 2 min)
#   Rscript 03-261012-gsea_l23.R --indir DIR --outdir DIR
#   Defaults: --indir and --outdir $M26_RESULTS_ROOT/cortex (or <repo>/results/cortex)
#
# Verification (2026-10-12)
#   Not reproducible exactly, because Reactome has grown since 2020: with reactome.db 1.89.0
#   and fgsea 1.32.2, 1,433 pathways are tested and 86 reach an adjusted P below 0.1 (2020:
#   8). The seven pathways of Fig. 3E that are still below 500 genes are among them, with
#   normalised enrichment scores within 0.1 of the 2020 values (1.42 to 1.64) and adjusted P
#   values of 0.010 to 0.11 (2020: 0.011 to 0.064); "Cell Cycle" now has more than 500 genes
#   and is not tested. The stored 2020 P values (down to 1.1e-4) imply about 10,000
#   permutations, not the 1,000 written in the archived code. Fig. 3E shows the stored 2020
#   result. With one value per gene instead of the most lengthened region, no pathway reaches
#   an adjusted P below 0.1.
#
# Provenance
#   Rewritten from human_asd_mapping.Rmd (2020 preprint of this work), which called
#   ReactomePA::gsePathway on the same ranking.
# =============================================================================

suppressMessages({
  library(dplyr)
  library(fgsea)
  library(reactome.db)
})

results_root <- Sys.getenv("M26_RESULTS_ROOT", "results")
opt <- list(indir = file.path(results_root, "cortex"), outdir = NULL)
args <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- args[i + 1]
}
if (is.null(opt$outdir)) opt$outdir <- opt$indir
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

load(file.path(opt$indir, "asd_utr_data.rda")) # mat_list, lm_tbl

m <- mat_list[["L2/3"]]
asd <- lm_tbl$Diagnose == "ASD"
score <- log(colMeans(m[asd, ], na.rm = TRUE) / colMeans(m[!asd, ], na.rm = TRUE))
names(score) <- sub("_.*", "", names(score)) # Entrez gene identifier of each 3'UTR region
score <- sort(score[is.finite(score)], decreasing = TRUE)
names(score) <- make.unique(names(score), sep = "_region") # lower-ranked regions of a gene: non-members

pathways <- as.list(reactomePATHID2EXTID)
pathways <- pathways[grepl("^R-HSA-", names(pathways))]
pathway_names <- unlist(as.list(reactomePATHID2NAME[names(pathways)]))
pathway_names <- sub("^Homo sapiens: ", "", pathway_names[!duplicated(names(pathway_names))])

set.seed(101)
res <- suppressWarnings(fgseaSimple(pathways, score, nperm = 10000, minSize = 10, maxSize = 500))
gsea <- as_tibble(res) |>
  transmute(
    ID = pathway, Description = unname(pathway_names[pathway]), setSize = size,
    enrichmentScore = ES, NES, pvalue = pval, p.adjust = padj
  ) |>
  arrange(p.adjust)

path_l23_tbl <- gsea |> filter(p.adjust < 0.1)
path_l23_tbl$Description <- factor(path_l23_tbl$Description, levels = rev(path_l23_tbl$Description[order(-path_l23_tbl$NES)]))
save(path_l23_tbl, file = file.path(opt$outdir, "asd_pathway.rda"))
write.table(gsea, file.path(opt$outdir, "gsea_l23_reactome.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(nrow(gsea), " pathways tested, ", nrow(path_l23_tbl), " with adjusted P < 0.1; written to ", opt$outdir)
