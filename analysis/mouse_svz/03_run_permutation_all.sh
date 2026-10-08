#!/usr/bin/env bash
# =============================================================================
# Driver for 260526-permutation_optimized_2.R (within-sample permutation LRT).
# Mouse SVZ.
#
# The mouse trajectory has a single lineage (no per-lineage split), so this
# is a one-shot dispatch.  Configured to match the existing 260401-...R LRT
# defaults (no subsampling -- mouse has ~14k lineage cells which is tractable
# at full size).
#
# Mouse-specific differences from the Knoblich driver:
#   - No --keep_induced_false_only (no `induced` column).
#   - One invocation (no lineage loop).
#   - SUBSAMPLE_N=0 by default (use full data; ~14k cells is fine at 26 cores).
#   - Pseudotime column is `dpt_pseudotime` (no lineage suffix).
#   - permutation_block_col is `sample_id` (6 samples; equivalent to the
#     "6 batches" mentioned in CLAUDE.md).
#
# Defaults are tuned for the mouse dataset on a 26-core box:
#   SUBSAMPLE_N=0        (no subsampling; mouse has ~14k cells which is fine)
#   NCORES=26
#   N_PERM_INIT=49 / N_PERM_MID=199 / N_PERM_MAX=999  (adaptive escalation)
#
# Outputs:
#   $OUTBASE/permutation_batch_number/
#     dirmult_regression_permutation_optimized2_*.csv
#     dirmult_regression_permutation_optimized2_*.rds
#     run_metadata_permutation_optimized2_*.json
#     cells_used.csv
#     run.log
#
# Usage:
#   bash run_permutation_all.sh
# or with overrides:
#   NCORES=16 SUBSAMPLE_N=8000 bash run_permutation_all.sh
# =============================================================================

set -euo pipefail
: "${M26_DATA_ROOT:?Set M26_DATA_ROOT to the data directory (see README)}"

# Pinned env
NCORES=${NCORES:-26}
SUBSAMPLE_N=${SUBSAMPLE_N:-0}
SEED=${SEED:-42}
SPLINE_DF=${SPLINE_DF:-5}
MIN_CELLS=${MIN_CELLS:-50}
N_PERM_INIT=${N_PERM_INIT:-49}
N_PERM_MID=${N_PERM_MID:-199}
N_PERM_MAX=${N_PERM_MAX:-999}
N_HITS_STOP=${N_HITS_STOP:-5}
PERMUTATION_STRATEGY=${PERMUTATION_STRATEGY:-within_sample}
PERMUTATION_BLOCK_COL=${PERMUTATION_BLOCK_COL:-sample_id}

# Inputs / outputs
INPUT=${INPUT:-${M26_DATA_ROOT}/3utr/wildtypes_analysis/2026_analyses/utrcalling/UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad}
OUTBASE=${OUTBASE:-${M26_DATA_ROOT}/3utr/wildtypes_analysis/2026_analyses/utrcalling/DirMult_optimized}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/03_260526-permutation_optimized_2.R"
RSCRIPT=${RSCRIPT:-Rscript}

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input h5ad not found: $INPUT" >&2
    exit 1
fi
if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: R script not found: $SCRIPT" >&2
    exit 1
fi

OUTDIR="$OUTBASE/permutation_batch_number"
mkdir -p "$OUTDIR"

echo "=========================================================================="
echo "Within-sample permutation LRT (260526) -- Mouse SVZ"
echo "  input:                 $INPUT"
echo "  outdir:                $OUTDIR"
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

"$RSCRIPT" "$SCRIPT" "$INPUT" \
    --pseudotime_col dpt_pseudotime \
    --covariate_cols batch_number \
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
echo "=========================================================================="
echo "Permutation run complete: $(date -Iseconds)"
echo "  Run the QC notebook to summarise concordance:"
echo "    $SCRIPT_DIR/260526-permutation_concordance.ipynb"
echo "=========================================================================="
