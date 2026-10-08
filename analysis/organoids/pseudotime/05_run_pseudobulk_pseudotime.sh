#!/usr/bin/env bash
# =============================================================================
# Driver for 260526-drimseq_pseudotime_pseudobulk.R.
#
# Pseudobulk DRIMSeq pseudotime LRT (controls-only) per lineage. Three
# quantile-bins of pseudotime as a factor; ~ super_batch + factor(pt_bin)
# vs ~ super_batch. Knoblich organoids.
#
# Outputs (per lineage):
#   $OUTBASE/<lineage>/pseudobulk_pseudotime/
#     drimseq_gene_<lineage>_pseudotime.csv
#     drimseq_feature_diagnostic_<lineage>_pseudotime.csv
#     delta_psi_pseudotime_gene_<lineage>.csv
#     delta_psi_pseudotime_feature_<lineage>.csv
#     drimseq_summary_<lineage>.csv
#     run_metadata_<lineage>.json
#     pseudobulk_units.csv
#     run.log
#
# Usage:
#   bash run_pseudobulk_pseudotime.sh
# or with overrides:
#   N_BINS=4 NCORES=8 bash run_pseudobulk_pseudotime.sh
# =============================================================================

set -euo pipefail
: "${M26_DATA_ROOT:?Set M26_DATA_ROOT to the data directory (see README)}"

# Pinned env
NCORES=${NCORES:-16}
N_BINS=${N_BINS:-3}
MIN_CELLS_PER_UNIT=${MIN_CELLS_PER_UNIT:-30}
MIN_UNITS_PER_BIN=${MIN_UNITS_PER_BIN:-4}
MIN_BIN_EXPR=${MIN_BIN_EXPR:-10}
MIN_SAMPS_GENE_EXPR=${MIN_SAMPS_GENE_EXPR:-5}
MIN_GENE_EXPR=${MIN_GENE_EXPR:-50}
MIN_SAMPS_FEATURE_EXPR=${MIN_SAMPS_FEATURE_EXPR:-5}
MIN_FEATURE_EXPR=${MIN_FEATURE_EXPR:-20}
SEED=${SEED:-42}
LINEAGES=${LINEAGES:-"inhibitory excitatory"}
# Batch covariate. Default: batch_number (= sample_id = "library" in Li et al.
# 2023; matches the cell-level LRT and permutation conventions). DO NOT change
# to super_batch unless you also re-read 260526-pseudobulk_pseudotime_methods.md
# Section 6 -- super_batch in controls-only has only 2 levels and one of them
# (213) contains a single sample, so it fails to control for between-sample
# variation among the other 7 controls samples.
COVARIATE_COL=${COVARIATE_COL:-batch_number}

# Inputs / outputs
INPUT=${INPUT:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad}
OUTBASE=${OUTBASE:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/dirmult_results}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/06_260526-drimseq_pseudotime_pseudobulk.R"
RSCRIPT=${RSCRIPT:-Rscript}

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input h5ad not found: $INPUT" >&2; exit 1
fi
if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: R script not found: $SCRIPT" >&2; exit 1
fi

echo "=========================================================================="
echo "Pseudobulk DRIMSeq pseudotime LRT (260526) -- Knoblich organoids"
echo "  input:                  $INPUT"
echo "  outbase:                $OUTBASE"
echo "  lineages:               $LINEAGES"
echo "  NCORES:                 $NCORES"
echo "  N_BINS:                 $N_BINS"
echo "  MIN_CELLS_PER_UNIT:     $MIN_CELLS_PER_UNIT"
echo "  MIN_UNITS_PER_BIN:      $MIN_UNITS_PER_BIN"
echo "  COVARIATE_COL:          $COVARIATE_COL  (matches Li et al. 2023 'library')"
echo "  SEED:                   $SEED"
echo "  start:                  $(date -Iseconds)"
echo "=========================================================================="

for lineage in $LINEAGES; do
    OUTDIR="$OUTBASE/$lineage/pseudobulk_pseudotime"
    mkdir -p "$OUTDIR"

    echo ""
    echo "--- Lineage: $lineage --- $(date -Iseconds) ---"
    echo "  output: $OUTDIR"
    echo ""

    "$RSCRIPT" "$SCRIPT" \
        --input "$INPUT" \
        --lineage "$lineage" \
        --pseudotime_col "dpt_pseudotime_${lineage}" \
        --covariate_col "$COVARIATE_COL" \
        --output_dir "$OUTDIR" \
        --keep_induced_false_only \
        --n_bins "$N_BINS" \
        --min_cells_per_unit "$MIN_CELLS_PER_UNIT" \
        --min_units_per_bin "$MIN_UNITS_PER_BIN" \
        --min_bin_expr "$MIN_BIN_EXPR" \
        --min_samps_gene_expr "$MIN_SAMPS_GENE_EXPR" \
        --min_gene_expr "$MIN_GENE_EXPR" \
        --min_samps_feature_expr "$MIN_SAMPS_FEATURE_EXPR" \
        --min_feature_expr "$MIN_FEATURE_EXPR" \
        --seed "$SEED" \
        --ncores "$NCORES" \
        2>&1 | tee "$OUTDIR/run.log"

    echo ""
    echo "--- Lineage $lineage finished: $(date -Iseconds) ---"
done

echo ""
echo "=========================================================================="
echo "All lineages complete: $(date -Iseconds)"
echo "  Run the QC notebook to compare against the cell-level LRT:"
echo "    $SCRIPT_DIR/260526-pseudobulk_concordance.ipynb"
echo "=========================================================================="
