#!/usr/bin/env bash
# =============================================================================
# Driver for 04_260526-permutation_optimized_2.R (within-sample permutation LRT).
#
# Loops over the two Knoblich lineages (inhibitory + excitatory) and runs the
# permutation LRT script on each. Same conventions as
# 01_run_dirmult_all.sh, with permutation-specific env vars added.
#
# Defaults are tuned for the Knoblich dataset on a 26-core box:
#   SUBSAMPLE_N=6000     (as in the LRT notebooks;
#                         0 = no subsampling)
#   NCORES=26
#   N_PERM_INIT=49 / N_PERM_MID=199 / N_PERM_MAX=999  (adaptive escalation)
#
# Outputs:
#   $OUTBASE/<lineage>/permutation_batch_number/
#     dirmult_regression_permutation_optimized2_*.csv
#     dirmult_regression_permutation_optimized2_*.rds
#     run_metadata_permutation_optimized2_*.json
#     cells_used.csv
#     run.log
#
# Usage:
#   bash 03_run_permutation_all.sh
# or with overrides:
#   NCORES=16 SUBSAMPLE_N=3000 bash 03_run_permutation_all.sh
# =============================================================================

set -euo pipefail
: "${M26_DATA_ROOT:?Set M26_DATA_ROOT to the data directory (see README)}"

# Pinned env (values used for the paper)
NCORES=${NCORES:-26}
SUBSAMPLE_N=${SUBSAMPLE_N:-6000}
SEED=${SEED:-42}
SPLINE_DF=${SPLINE_DF:-5}
MIN_CELLS=${MIN_CELLS:-50}
N_PERM_INIT=${N_PERM_INIT:-49}
N_PERM_MID=${N_PERM_MID:-199}
N_PERM_MAX=${N_PERM_MAX:-999}
N_HITS_STOP=${N_HITS_STOP:-5}
PERMUTATION_STRATEGY=${PERMUTATION_STRATEGY:-within_sample}
PERMUTATION_BLOCK_COL=${PERMUTATION_BLOCK_COL:-sample_id}
LINEAGES=${LINEAGES:-"inhibitory excitatory"}

# Inputs / outputs
INPUT=${INPUT:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad}
OUTBASE=${OUTBASE:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/dirmult_results}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/04_260526-permutation_optimized_2.R"
RSCRIPT=${RSCRIPT:-Rscript}

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input h5ad not found: $INPUT" >&2
    exit 1
fi
if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: R script not found: $SCRIPT" >&2
    exit 1
fi

echo "=========================================================================="
echo "Within-sample permutation LRT (260526) -- Knoblich organoids"
echo "  input:                 $INPUT"
echo "  outbase:               $OUTBASE"
echo "  lineages:              $LINEAGES"
echo "  NCORES:                $NCORES"
echo "  SUBSAMPLE_N:           $SUBSAMPLE_N"
echo "  SEED:                  $SEED"
echo "  SPLINE_DF:             $SPLINE_DF"
echo "  MIN_CELLS:             $MIN_CELLS"
echo "  N_PERM_INIT / MID / MAX: $N_PERM_INIT / $N_PERM_MID / $N_PERM_MAX"
echo "  N_HITS_STOP:           $N_HITS_STOP"
echo "  PERMUTATION_STRATEGY:  $PERMUTATION_STRATEGY"
echo "  PERMUTATION_BLOCK_COL: $PERMUTATION_BLOCK_COL"
echo "  start:                 $(date -Iseconds)"
echo "=========================================================================="

for lineage in $LINEAGES; do
    OUTDIR="$OUTBASE/$lineage/permutation_batch_number"
    mkdir -p "$OUTDIR"

    echo ""
    echo "--- Lineage: $lineage --- $(date -Iseconds) ---"
    echo "  output: $OUTDIR"
    echo ""

    "$RSCRIPT" "$SCRIPT" "$INPUT" \
        --pseudotime_col "dpt_pseudotime_${lineage}" \
        --covariate_cols batch_number \
        --keep_induced_false_only \
        --permutation_block_col "$PERMUTATION_BLOCK_COL" \
        --permutation_strategy "$PERMUTATION_STRATEGY" \
        --subsample_n "$SUBSAMPLE_N" \
        --seed "$SEED" \
        --output_dir "$OUTDIR" \
        --ncores "$NCORES" \
        --spline_df "$SPLINE_DF" \
        --min_cells "$MIN_CELLS" \
        --n_perm_init "$N_PERM_INIT" \
        --n_perm_mid "$N_PERM_MID" \
        --n_perm_max "$N_PERM_MAX" \
        --n_hits_stop "$N_HITS_STOP" \
        --calibration_n 200 \
        2>&1 | tee "$OUTDIR/run.log"

    echo ""
    echo "--- Lineage $lineage finished: $(date -Iseconds) ---"
done

echo ""
echo "=========================================================================="
echo "All lineages complete: $(date -Iseconds)"
echo "  Run the QC notebook to summarise concordance:"
echo "    $SCRIPT_DIR/260526-permutation_concordance.ipynb"
echo "=========================================================================="
