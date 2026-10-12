#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# 01-261012-utr_lengths_per_cell_type.R
#
# Purpose
#   3'UTR lengths of the post-mortem cortex samples (Fig. 3), from the per-molecule
#   positions of 03-261012-read_3utr_positions.R and the cell-type labels of
#   Velmeshev et al. 2019:
#   - samples: 10x libraries of ASD and control individuals younger than 30 years
#     (23 prefrontal and 18 anterior cingulate samples from 31 individuals);
#   - per gene: for each sample and cell type, the mean length over the nuclei of that
#     type, counting one molecule per nucleus (the first in barcode-UMI order); genes
#     are kept if they have positions in every sample;
#   - per cell type: for each sample, the mean length over all 3'UTR molecules of the
#     nuclei of that type (Fig. 3B).
#
# Inputs (<cortex_dir>)
#   prefrontal_cortex.rda, anterior_cortex.rda   per-molecule positions (prefrontal_list, anterior_list)
#   srr_annotation.csv                           SRA run table of BioProject PRJNA434002
#                                                (semicolon-separated, no header)
#   cell_annotation.csv                          per-nucleus labels of Velmeshev et al. 2019
#                                                (semicolon-separated, decimal comma)
#
# Outputs (<outdir>)
#   utr_per_cell_type_matrix.rda   mat_list (one samples x genes matrix per cell type), lm_tbl (samples)
#   meta_utr_asd_.rda              meta_utr_asd_tbl (one row per sample and cell type; Cell_Type
#                                  levels ordered by mean length, as plotted in Fig. 3B)
#
# How to run (conda env apa-r; about 20 min, 12 GB of memory)
#   Rscript 01-261012-utr_lengths_per_cell_type.R --cortex_dir DIR --outdir DIR
#   Defaults: --cortex_dir $M26_DATA_ROOT/copied_manuel_data/Mapping_3UTR_10X,
#             --outdir $M26_RESULTS_ROOT/cortex (or <repo>/results/cortex)
#
# Verification (2026-10-12)
#   Run on the deposited position lists, mat_list (all 17 cell types, values, missing
#   values and dimnames) and lm_tbl (Diagnose, ID, Region, Sex, Age) are identical to the
#   stored asd_utr_data.rda, and every value of meta_utr_asd_tbl is identical to the stored
#   meta_utr_asd_.rda (the rows are in a different order). 21 min, 12 GB of memory.
#
# Provenance
#   Rewritten from the 2019-2021 R code of the 2020 preprint (human_asd_mapping.Rmd for the
#   per-gene matrices; NEW_ASD_Analysis.Rmd for the per-cell-type means).
# =============================================================================

suppressMessages({
  library(dplyr)
  library(readr)
})

data_root <- Sys.getenv("M26_DATA_ROOT")
results_root <- Sys.getenv("M26_RESULTS_ROOT", "results")
opt <- list(
  cortex_dir = file.path(data_root, "copied_manuel_data", "Mapping_3UTR_10X"),
  outdir = file.path(results_root, "cortex")
)
args <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opt)) stop("unknown argument: ", args[i])
  opt[[key]] <- args[i + 1]
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- samples -----------------------------------------------------------------
srr <- read_csv2(file.path(opt$cortex_dir, "srr_annotation.csv"), col_names = FALSE, show_col_types = FALSE)
runs <- srr |>
  filter(grepl("10x", X12), X15 %in% c("ASD", "Control")) |>
  transmute(ID = substring(X12, 1, 4), Region = X17, Age = as.integer(X14), Sex = X18, Diagnose = X15)
pfc <- runs |>
  filter(Region == "Brain, prefrontal cortex") |>
  distinct(ID, .keep_all = TRUE)
acc <- runs |> filter(Region == "Brain, anterior cingulate cortex")
lm_tbl <- bind_rows(pfc, acc) |> filter(Age < 30)
lm_tbl$area <- ifelse(lm_tbl$Region == "Brain, prefrontal cortex", "PFC", "ACC")
message(nrow(lm_tbl), " samples from ", n_distinct(lm_tbl$ID), " individuals")

# ---- positions and cell-type labels --------------------------------------------
load(file.path(opt$cortex_dir, "prefrontal_cortex.rda"))
load(file.path(opt$cortex_dir, "anterior_cortex.rda"))
names(prefrontal_list) <- sub(".*RUN_", "", names(prefrontal_list))
names(anterior_list) <- sub(".*RUN_", "", names(anterior_list))
positions <- c(
  prefrontal_list[lm_tbl$ID[lm_tbl$area == "PFC"]],
  anterior_list[lm_tbl$ID[lm_tbl$area == "ACC"]]
)
rm(prefrontal_list, anterior_list)
stopifnot(!any(vapply(positions, is.null, logical(1))))

anno <- read_csv2(file.path(opt$cortex_dir, "cell_annotation.csv"), show_col_types = FALSE) |>
  filter(!is.na(cluster)) |>
  transmute(individual = as.character(individual), area = region, barcode = substring(cell, 1, 18), cluster)
# cell types in order of first appearance in the prefrontal samples (the order of mat_list)
pfc_ids <- lm_tbl$ID[lm_tbl$area == "PFC"]
cell_types <- unique(unlist(lapply(pfc_ids, function(id) anno$cluster[anno$individual == id & anno$area == "PFC"])))

# per sample: mean length per gene and cell type (one molecule per nucleus), and the
# mean over all molecules per cell type
per_sample <- function(gene_positions, id, area) {
  lab <- anno[anno$individual == id & anno$area == area, ]
  pos <- unlist(gene_positions, use.names = FALSE)
  barcode <- unlist(lapply(gene_positions, names), use.names = FALSE)
  gene <- rep(names(gene_positions), lengths(gene_positions))
  cluster <- lab$cluster[match(barcode, lab$barcode)]

  meta <- tapply(pos[!is.na(cluster)], factor(cluster[!is.na(cluster)], levels = cell_types), mean)

  first <- !duplicated(paste(gene, barcode)) & !is.na(cluster)
  key <- factor(paste(gene[first], cluster[first], sep = "\t"))
  m <- tapply(pos[first], key, mean)
  parts <- do.call(rbind, strsplit(names(m), "\t", fixed = TRUE))
  gene_means <- matrix(NA_real_, length(gene_positions), length(cell_types),
    dimnames = list(names(gene_positions), cell_types)
  )
  gene_means[cbind(match(parts[, 1], names(gene_positions)), match(parts[, 2], cell_types))] <- m
  list(meta = meta, gene_means = gene_means)
}

res <- Map(per_sample, positions, lm_tbl$ID, lm_tbl$area)
genes <- Reduce(intersect, lapply(positions, names))
mat_list <- lapply(setNames(cell_types, cell_types), function(ct) {
  m <- t(vapply(res, function(r) r$gene_means[genes, ct], numeric(length(genes))))
  dimnames(m) <- list(lm_tbl$ID, genes)
  m
})

meta_utr_asd_tbl <- bind_rows(lapply(seq_along(res), function(i) {
  tibble(
    Diagnose = lm_tbl$Diagnose[i], Individual = lm_tbl$ID[i], Cell_Type = cell_types,
    UTR3_Length = as.numeric(res[[i]]$meta[cell_types]), Region = lm_tbl$area[i]
  )
}))
level_order <- meta_utr_asd_tbl |>
  group_by(Cell_Type) |>
  summarise(m = mean(UTR3_Length, na.rm = TRUE)) |>
  arrange(m)
meta_utr_asd_tbl$Cell_Type <- factor(meta_utr_asd_tbl$Cell_Type, levels = level_order$Cell_Type)
meta_utr_asd_tbl$Diagnose <- factor(meta_utr_asd_tbl$Diagnose)

lm_tbl <- lm_tbl |> select(Diagnose, ID, Region, Age, Sex)
save(mat_list, lm_tbl, file = file.path(opt$outdir, "utr_per_cell_type_matrix.rda"))
save(meta_utr_asd_tbl, file = file.path(opt$outdir, "meta_utr_asd_.rda"))
message(length(genes), " genes, ", length(cell_types), " cell types; written to ", opt$outdir)
