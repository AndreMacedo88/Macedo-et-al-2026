#!/bin/bash
# ============================================================================
# Run Dirichlet-multinomial (optimized_2) regression on Knoblich CONTROLS ONLY
# (induced=False), one per lineage, with batch_number as fixed effect.
#
# Uses 260503-dirichlet_multinomial_analysis_optimized_2.R (publication-rigor:
# shared-s, conditional-MLE refits, padded-null safety net, optional
# parametric bootstrap, real-data null calibration, 16-test verification).
#
# Rationale for the analysis design:
#   * `induced` is perfectly confounded with `sample`/`batch_number` in the
#     Knoblich data (every sample is 100% one induced state). Restricting to
#     induced=False removes the confound, so batch_number becomes a clean
#     8-level fixed effect.
#   * Pseudotime was computed without batch correction; including batch_number
#     in the regression now correctly absorbs technical batch effects.
#   * The full inhibitory induced=False pool has 5,130 cells, excitatory
#     16,153. Both subsampled to a fixed N for tractability and
#     reproducibility (seed-controlled). cells_used.csv is written into each
#     run folder for full reproducibility.
#
# Usage:
#   bash run_dirmult_all.sh                          # defaults (N=2000)
#   NCORES=26 SUBSAMPLE_N=2500 bash run_dirmult_all.sh
#   NCORES=26 BOOTSTRAP=199 bash run_dirmult_all.sh  # enable bootstrap
#
# Tunables (env vars):
#   NCORES          cores for parallel per-gene fitting     (default 16)
#   SUBSAMPLE_N     cells kept per lineage after filters    (default 2000)
#   SEED            RNG seed for the subsample              (default 42)
#   SPLINE_DF       df for bs() spline of pseudotime        (default 5)
#   MIN_CELLS       skip UTRs with fewer than this many     (default 50)
#   BOOTSTRAP       N bootstrap reps; 0 = off               (default 0)
#   BOOTSTRAP_PVAL  bootstrap genes with p_chisq < this     (default 0.1)
#   CALIBRATION_N   real-data null calibration sample       (default 200)
#   INPUT           path to the raw-UTR + pseudotime h5ad
#   RESULTS_ROOT    root output directory
#
# Outputs under <RESULTS_ROOT>/<lineage>/controls_only_batch_number/
#   dirmult_regression_*.csv  per-gene LRT statistics + padj (with shared-s,
#                             ll_full/ll_null, lrt_status, max/min_var_ratio,
#                             cells_per_param, full_converged columns)
#   dirmult_regression_*.rds  full output (results, SE, trend, bootstrap,
#                             calibration, opts)
#   cells_used.csv            cell barcodes + pseudotime + induced + batch
#   run.log                   stdout/stderr of the run
#
# Runtime estimates (method_2 is ~7x slower per gene than method_1 because of
# the conditional-MLE refit and variance-ratio diagnostics).
# Empirical: 200 genes / N=1000 cells / 4 cores = 7 min wallclock.
# Per-gene scales as roughly N^1.72.
# Full data = 14,755 UTRs.
#   N=1500 cells, 26 cores -> ~3h per lineage, ~6h  total ✓ short overnight
#   N=2000 cells, 26 cores -> ~4.5h/lineage, ~9h   total ✓ overnight (default)
#   N=2500 cells, 26 cores -> ~6.5h/lineage, ~13h  total (long overnight)
#   N=3000 cells, 26 cores -> ~9h /lineage,  ~18h  total (>1 night)
# Phase D (real-data null calibration, default n=200) adds ~30s on the same
# core count. Bootstrap (off by default) is the biggest opt-in cost: with
# B=199 on ~5-10% borderline genes, expect a 30-60% wallclock increase.
# ============================================================================

set -euo pipefail
: "${M26_DATA_ROOT:?Set M26_DATA_ROOT to the data directory (see README)}"

# --- Configuration (override via environment variables) ---------------------
INPUT="${INPUT:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad}"
RSCRIPT="${RSCRIPT:-Rscript}"
SCRIPT="${SCRIPT:-$(dirname "$(readlink -f "$0")")/02_260503-dirichlet_multinomial_analysis_optimized_2.R}"
RESULTS_ROOT="${RESULTS_ROOT:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/dirmult_results}"
NCORES="${NCORES:-26}"
SUBSAMPLE_N="${SUBSAMPLE_N:-6000}"   # value used for the paper
SEED="${SEED:-42}"
SPLINE_DF="${SPLINE_DF:-5}"
MIN_CELLS="${MIN_CELLS:-50}"
BOOTSTRAP="${BOOTSTRAP:-199}"          # value used for the paper; 0=off
BOOTSTRAP_PVAL="${BOOTSTRAP_PVAL:-0.1}" # bootstrap genes with p_chisq < this
CALIBRATION_N="${CALIBRATION_N:-200}"   # real-data null calibration footprint

# --- Sanity checks ----------------------------------------------------------
[ -f "$INPUT" ]   || { echo "Input not found: $INPUT" >&2; exit 1; }
command -v "$RSCRIPT" >/dev/null || { echo "Rscript not found or not executable: $RSCRIPT" >&2; exit 1; }
[ -f "$SCRIPT" ]  || { echo "R script not found: $SCRIPT" >&2; exit 1; }

mkdir -p "$RESULTS_ROOT"

echo "============================================================"
echo "Dirichlet-multinomial regression (optimized_2): Knoblich controls only"
echo "  input:          $INPUT"
echo "  script:         $SCRIPT"
echo "  results:        $RESULTS_ROOT"
echo "  ncores:         $NCORES"
echo "  subsample_n:    $SUBSAMPLE_N"
echo "  seed:           $SEED"
echo "  spline_df:      $SPLINE_DF"
echo "  min_cells:      $MIN_CELLS"
echo "  bootstrap:      $BOOTSTRAP"
echo "  bootstrap_pval: $BOOTSTRAP_PVAL"
echo "  calibration_n:  $CALIBRATION_N"
echo "  start:          $(date)"
echo "============================================================"

run_one() {
    local lineage="$1"
    local outdir="$RESULTS_ROOT/$lineage/controls_only_batch_number"

    mkdir -p "$outdir"
    local logfile="$outdir/run.log"

    echo ""
    echo "------------------------------------------------------------"
    echo "[$(date)] $lineage / controls_only_batch_number"
    echo "  output: $outdir"
    echo "  log:    $logfile"
    echo "------------------------------------------------------------"

    "$RSCRIPT" "$SCRIPT" "$INPUT" \
        --pseudotime_col "dpt_pseudotime_${lineage}" \
        --covariate_cols batch_number \
        --keep_induced_false_only \
        --subsample_n "$SUBSAMPLE_N" \
        --seed "$SEED" \
        --output_dir "$outdir" \
        --ncores "$NCORES" \
        --spline_df "$SPLINE_DF" \
        --min_cells "$MIN_CELLS" \
        --bootstrap "$BOOTSTRAP" \
        --bootstrap_pval "$BOOTSTRAP_PVAL" \
        --calibration_n "$CALIBRATION_N" \
        2>&1 | tee "$logfile"
}

# Inhibitory first (smaller cohort if not subsampled; here both end up at
# SUBSAMPLE_N so order is just a convention).
for lineage in inhibitory excitatory; do
    run_one "$lineage"
done

echo ""
echo "============================================================"
echo "[$(date)] Both runs complete."
echo "Results under: $RESULTS_ROOT"
echo "============================================================"
ls -la "$RESULTS_ROOT"/*/controls_only_batch_number/dirmult_regression_*.csv 2>/dev/null || true
