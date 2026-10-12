#!/usr/bin/env -S Rscript --vanilla
# =============================================================================
# 01-261012-utr_reference.R
#
# Purpose
#   3'UTR regions of the human genome for the post-mortem cortex data (Fig. 3).
#   The 3'UTRs of all transcripts in TxDb.Hsapiens.UCSC.hg38.knownGene are merged
#   into non-overlapping regions, each region is named after the Entrez gene it
#   overlaps (a gene with several separate regions gets the suffixes _1, _2, ...),
#   and only autosomal regions longer than 100 bp are kept. Chromosome names lose
#   the "chr" prefix to match the Cell Ranger GRCh38 reference.
#
# Inputs
#   TxDb.Hsapiens.UCSC.hg38.knownGene (Bioconductor annotation package)
#
# Outputs
#   <outdir>/threeUTR_human.rda   GRanges `threeUTR_human`
#
# How to run (conda env apa-r with GenomicFeatures and the TxDb package)
#   Rscript 01-261012-utr_reference.R --outdir /path/to/cortex
#
# Provenance
#   Written in 2019 for the 2020 preprint of this work; the version of the TxDb
#   package used then was not recorded, so a current version can give slightly
#   different regions. The regions used for the paper are deposited with the
#   processed cortex data (threeUTR_human.rda).
# =============================================================================

suppressMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) == 2 && args[1] == "--outdir") args[2] else stop("usage: --outdir DIR")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene <- transcriptsBy(txdb, "gene")
utrs <- unlist(as(threeUTRsByTranscript(txdb), "GRangesList"))
regions <- reduce(utrs)

# name each merged region after the gene it overlaps (the last overlapping gene wins)
hits <- findOverlaps(regions, tx_by_gene)
names(regions)[queryHits(hits)] <- names(tx_by_gene)[subjectHits(hits)]
regions <- regions[!is.na(names(regions))]
names(regions) <- make.unique(names(regions), sep = "_")

seqlevels(regions) <- sub("^chr", "", seqlevels(regions))
regions <- keepSeqlevels(regions, as.character(1:22), pruning.mode = "coarse")
threeUTR_human <- regions[end(regions) - start(regions) > 100]

save(threeUTR_human, file = file.path(outdir, "threeUTR_human.rda"))
message(length(threeUTR_human), " regions written to ", outdir)
