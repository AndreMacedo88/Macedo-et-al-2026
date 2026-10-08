#!/bin/bash
# =============================================================================
# Run the primary treatment-effect screen (main-only, no interaction term).
# One invocation per lineage (inhibitory, excitatory).
#
# Tests one LRT per gene per lineage:
#   full ~ bs(pt, df) + super_batch + induced
#   null ~ bs(pt, df) + super_batch
# Reports the `induced` coefficient as the treatment-effect endpoint
# (constant treatment shift in UTR usage at any pseudotime point).
#
# Per-gene precision: gene-wise (DRIMSeq default) with NaN/Inf/<=0 entries
# repaired to common precision (typically affects ~150 high-K genes out of
# ~5,000-7,500 post-filter).
#
# Settings (env vars): see the configuration block below. Defaults are the
# values used for the paper.
#
# Usage:
#   bash 01_run_treatment_main.sh                         # defaults
#   NCORES=26 SPLINE_DF=5 N_BINS=8 bash 01_run_treatment_main.sh
#   BOOTSTRAP_B=199 NCORES=26 bash 01_run_treatment_main.sh   # supplementary
#
# The induced arm is restricted to cells with a targeting guide
# (--require_targeting_guide TRUE); uninduced cells are kept regardless of guide.
#
# Outputs under <RESULTS_ROOT>/<lineage_slug>/main_only/:
#   drimseq_gene_<lineage>_treatment_main.csv         (primary, single LRT)
#   drimseq_feature_diagnostic_<lineage>_treatment_main.csv
#   delta_psi_treatment_gene_<lineage>.csv
#   delta_psi_treatment_feature_<lineage>.csv
#   pseudobulk_units.csv
#   drimseq_summary_<lineage>.csv
#   bootstrap_stability_gene_<lineage>.csv      (if BOOTSTRAP_B > 0)
#   bootstrap_stability_summary_<lineage>.csv   (if BOOTSTRAP_B > 0)
#   run.log
#
# Runtime estimates:
#   Per-lineage observed run (BOOTSTRAP_B=0):       ~10-25 min
#   With bootstrap (BOOTSTRAP_B=49):                ~6-10 hours per lineage
#   With bootstrap (BOOTSTRAP_B=199):               ~24-30 hours per lineage
# =============================================================================

set -euo pipefail
: "${M26_DATA_ROOT:?Set M26_DATA_ROOT to the data directory (see README)}"

# --- Configuration -----------------------------------------------------------
INPUT="${INPUT:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/merged_bw12_above10cts_apa_named_annotated_fix_pt.h5ad}"
RSCRIPT="${RSCRIPT:-Rscript}"
SCRIPT="${SCRIPT:-$(dirname "$(readlink -f "$0")")/02_drimseq_treatment_main_pseudobulk.R}"
RESULTS_ROOT="${RESULTS_ROOT:-${M26_DATA_ROOT}/3utr/organoids/knoblich/utrcalling/dirmult_treatment_results}"
NCORES="${NCORES:-16}"
SPLINE_DF="${SPLINE_DF:-5}"
N_BINS="${N_BINS:-8}"
MIN_CELLS_PER_UNIT="${MIN_CELLS_PER_UNIT:-10}"
MIN_SAMPS_GENE_EXPR="${MIN_SAMPS_GENE_EXPR:-5}"
MIN_GENE_EXPR="${MIN_GENE_EXPR:-50}"
MIN_SAMPS_FEATURE_EXPR="${MIN_SAMPS_FEATURE_EXPR:-5}"
MIN_FEATURE_EXPR="${MIN_FEATURE_EXPR:-20}"
MIN_UNITS_PER_COND="${MIN_UNITS_PER_COND:-3}"
MIN_COND_EXPR="${MIN_COND_EXPR:-10}"
BOOTSTRAP_B="${BOOTSTRAP_B:-0}"
BOOTSTRAP_SEED="${BOOTSTRAP_SEED:-2026}"
SEED="${SEED:-42}"

[ -f "$INPUT" ]   || { echo "Input not found: $INPUT" >&2; exit 1; }
command -v "$RSCRIPT" >/dev/null || { echo "Rscript not found: $RSCRIPT" >&2; exit 1; }
[ -f "$SCRIPT" ]  || { echo "R script not found: $SCRIPT" >&2; exit 1; }
mkdir -p "$RESULTS_ROOT"

echo "============================================================"
echo "DRIMSeq treatment MAIN-EFFECT screen (no interaction)"
echo "  input:                $INPUT"
echo "  script:               $SCRIPT"
echo "  results:              $RESULTS_ROOT"
echo "  ncores:               $NCORES   spline_df: $SPLINE_DF   n_bins: $N_BINS"
echo "  min_cells_per_unit:   $MIN_CELLS_PER_UNIT"
echo "  min_units_per_cond:   $MIN_UNITS_PER_COND   min_cond_expr: $MIN_COND_EXPR"
echo "  bootstrap_b:          $BOOTSTRAP_B   seed: $SEED"
echo "  start:                $(date)"
echo "============================================================"

run_one() {
    local lineage="$1"
    local slug
    slug="$(echo "$lineage" | tr '[:upper:] ' '[:lower:]_')"
    local outdir="$RESULTS_ROOT/$slug/main_only"
    mkdir -p "$outdir"
    local logfile="$outdir/run.log"

    echo ""
    echo "------------------------------------------------------------"
    echo "[$(date)] lineage: $lineage  (slug: $slug)"
    echo "  output: $outdir"
    echo "  log:    $logfile"
    echo "------------------------------------------------------------"

    "$RSCRIPT" "$SCRIPT" \
        --input "$INPUT" \
        --lineage "$lineage" \
        --output_dir "$outdir" \
        --spline_df "$SPLINE_DF" \
        --n_bins "$N_BINS" \
        --min_cells_per_unit "$MIN_CELLS_PER_UNIT" \
        --min_samps_gene_expr "$MIN_SAMPS_GENE_EXPR" \
        --min_gene_expr "$MIN_GENE_EXPR" \
        --min_samps_feature_expr "$MIN_SAMPS_FEATURE_EXPR" \
        --min_feature_expr "$MIN_FEATURE_EXPR" \
        --min_units_per_cond "$MIN_UNITS_PER_COND" \
        --min_cond_expr "$MIN_COND_EXPR" \
        --bootstrap_b "$BOOTSTRAP_B" \
        --bootstrap_seed "$BOOTSTRAP_SEED" \
        --seed "$SEED" \
        --ncores "$NCORES" \
        --require_targeting_guide TRUE \
        2>&1 | tee "$logfile"
}

for lineage in inhibitory excitatory; do
    run_one "$lineage"
done

echo ""
echo "============================================================"
echo "[$(date)] Both lineages complete (main-effect screen)."
echo "Results under: $RESULTS_ROOT/<lineage>/main_only/"
echo "============================================================"
ls -la "$RESULTS_ROOT"/*/main_only/drimseq_gene_*_treatment_main.csv 2>/dev/null || true
