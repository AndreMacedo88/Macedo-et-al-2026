#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# 03-261012-read_3utr_positions.R
#
# Purpose
#   3'-end position of every molecule in each 3'UTR region, per nucleus, for the
#   Cell Ranger runs of one brain region (Fig. 3). For each alignment in a region
#   (mapping quality >= 10, no duplicates or secondary alignments) the position of
#   its 3' end is measured from the start of the region on the gene's strand; reads
#   that share a cell barcode and UMI are reduced to their median position. Positions
#   more than 100 nt past the end of the region (spliced or read-through alignments)
#   are dropped, and regions with a single molecule in the run are skipped.
#
# Inputs
#   <runs_dir>/RUN_<individual>/outs/possorted_genome_bam.bam   (02-cellranger_count.sh)
#   <utr_rda>   threeUTR_human.rda                              (01-261012-utr_reference.R)
#
# Outputs
#   <out_rda>   a list with one element per run, named "<label>_<run directory>";
#               each element is a list with one element per 3'UTR region, a numeric
#               vector of positions (bp) named by cell barcode, one value per molecule.
#               The paper used label prefrontal_cortex (prefrontal_cortex.rda, object
#               prefrontal_list) and anterior_cortex (anterior_cortex.rda, object anterior_list).
#
# How to run (conda env apa-r with GenomicAlignments, Rsamtools and BiocParallel)
#   Rscript 03-261012-read_3utr_positions.R RUNS_DIR UTR_RDA LABEL OUT_RDA [WORKERS]
#   e.g. Rscript 03-261012-read_3utr_positions.R pfc_runs threeUTR_human.rda prefrontal_cortex \
#        prefrontal_cortex.rda 16
#
# Provenance
#   Written in 2019 for the 2020 preprint of this work (read_in_10Xutr and
#   extract3utrPositions10X); unused columns (sequence, CIGAR, read length) were removed.
#   The BAM files are no longer available, so this script could not be re-run; the
#   position lists used for the paper are deposited with the processed cortex data.
# =============================================================================

suppressMessages({
  library(BiocParallel)
  library(GenomicAlignments)
  library(Rsamtools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: RUNS_DIR UTR_RDA LABEL OUT_RDA [WORKERS]")
runs_dir <- args[1]
load(args[2]) # threeUTR_human
label <- args[3]
out_rda <- args[4]
workers <- if (length(args) >= 5) as.integer(args[5]) else 16L
chunk_size <- 250

# positions of the molecules in one chunk of regions of one BAM file
chunk_positions <- function(idx, bam, regions) {
  param <- ScanBamParam(
    flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
    mapqFilter = 10, which = regions[idx], tag = c("CB", "UB")
  )
  aln <- readGAlignments(bam, param = param, with.which_label = TRUE)
  if (length(aln) == 0) {
    return(list())
  }
  region_of <- as.character(mcols(aln)$which_label)
  by_region <- split(seq_along(aln), factor(region_of, levels = unique(region_of)))
  labels <- sprintf("%s:%d-%d", seqnames(regions[idx]), start(regions[idx]), end(regions[idx]))
  out <- list()
  for (lab in names(by_region)) {
    r <- regions[idx][match(lab, labels)]
    a <- aln[by_region[[lab]]]
    pos <- if (as.character(strand(r)) == "+") end(a) - start(r) else end(r) - start(a)
    molecule <- paste0(mcols(a)$CB, "_", mcols(a)$UB)
    per_molecule <- tapply(pos, molecule, median)
    names(per_molecule) <- substring(names(per_molecule), 1, 18) # cell barcode, e.g. AAACCTGAGAAACCAT-1
    per_molecule <- per_molecule[which(per_molecule < width(r) + 100)]
    if (length(per_molecule) > 1) out[[names(r)]] <- per_molecule
  }
  out
}

run_positions <- function(run) {
  bam <- file.path(run, "outs", "possorted_genome_bam.bam")
  chunks <- split(seq_along(threeUTR_human), ceiling(seq_along(threeUTR_human) / chunk_size))
  res <- bplapply(chunks, chunk_positions,
    bam = bam, regions = threeUTR_human,
    BPPARAM = SnowParam(workers = workers, type = "SOCK")
  )
  unlist(unname(res), recursive = FALSE)
}

runs <- list.files(runs_dir, pattern = "^RUN_", full.names = TRUE)
runs <- runs[file.exists(file.path(runs, "outs", "possorted_genome_bam.bam"))]
positions <- lapply(runs, run_positions)
names(positions) <- paste0(label, "_", runs)

obj <- paste0(sub("_cortex$", "", label), "_list")
assign(obj, positions)
save(list = obj, file = out_rda)
message(length(positions), " runs written to ", out_rda, " (object ", obj, ")")
