#!/usr/bin/env bash
# Cell Ranger 2.2.0 alignment of the post-mortem cortex libraries (Velmeshev et al. 2019,
# Sequence Read Archive BioProject PRJNA434002), one run per library.
#
# Usage: 02-cellranger_count.sh FASTQ_DIR OUT_DIR SAMPLE [SAMPLE ...]
#   FASTQ_DIR  folder with the FASTQ files as submitted to the SRA (e.g. 5864_BA9_S14_L001_R2_001.fastq.gz)
#   OUT_DIR    folder in which the runs RUN_<individual> are written
#   SAMPLE     FASTQ sample name, e.g. 5864_BA9 (prefrontal) or 5864_BA24 (anterior cingulate)
# The reference is set with CELLRANGER_REF (default refdata-cellranger-GRCh38-3.0.0).
# Prefrontal and anterior cingulate libraries go to separate OUT_DIRs, because the
# runs are named by individual only.
set -euo pipefail

fastq_dir=$1
out_dir=$2
shift 2
ref=${CELLRANGER_REF:-refdata-cellranger-GRCh38-3.0.0}

mkdir -p "$out_dir"
cd "$out_dir"
for sample in "$@"; do
    cellranger count --id="RUN_${sample%%_*}" --sample="$sample" --fastqs="$fastq_dir" --transcriptome="$ref"
done
