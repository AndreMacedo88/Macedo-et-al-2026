#!/usr/bin/env bash
# =============================================================================
# Driver for 04_260526-drimseq_pseudotime_pseudobulk.R.
# Mouse SVZ. Single lineage (no per-lineage loop, no controls filter).
#
# Outputs:
#   $OUTBASE/pseudobulk_pseudotime/
#     drimseq_gene_all_pseudotime.csv
#     drimseq_feature_diagnostic_all_pseudotime.csv
#     delta_psi_pseudotime_gene_all.csv
#     delta_psi_pseudotime_feature_all.csv
#     drimseq_summary_all.csv
#     run_metadata_all.json
#     pseudobulk_units.csv
#     run.log
#
# Usage:
#   bash 04_run_pseudobulk_pseudotime.sh
# or with overrides:
#   N_BINS=4 NCORES=8 bash 04_run_pseudobulk_pseudotime.sh
# =============================================================================

set -euo pipefail
: "${M26_DATA_ROOT:?Set M26_DATA_ROOT to the data directory (see README)}"

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
PSEUDOTIME_COL=${PSEUDOTIME_COL:-dpt_pseudotime}
COVARIATE_COL=${COVARIATE_COL:-batch_number}

# Inputs / outputs
INPUT=${INPUT:-${M26_DATA_ROOT}/3utr/wildtypes_analysis/2026_analyses/utrcalling/UTR_count_table_genelevel_bw12_above10cts_apa_named.h5ad}
OUTBASE=${OUTBASE:-${M26_DATA_ROOT}/3utr/wildtypes_analysis/2026_analyses/utrcalling/DirMult_optimized}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/04_260526-drimseq_pseudotime_pseudobulk.R"
RSCRIPT=${RSCRIPT:-Rscript}

if [[ ! -f "$INPUT" ]]; then echo "ERROR: input h5ad not found: $INPUT" >&2; exit 1; fi
if [[ ! -f "$SCRIPT" ]]; then echo "ERROR: R script not found: $SCRIPT" >&2; exit 1; fi

OUTDIR="$OUTBASE/pseudobulk_pseudotime"
mkdir -p "$OUTDIR"

echo "=========================================================================="
echo "Pseudobulk DRIMSeq pseudotime LRT (260526) -- Mouse SVZ"
echo "  input:                  $INPUT"
echo "  outdir:                 $OUTDIR"
echo "  PSEUDOTIME_COL:         $PSEUDOTIME_COL"
echo "  COVARIATE_COL:          $COVARIATE_COL"
echo "  NCORES:                 $NCORES"
echo "  N_BINS:                 $N_BINS"
echo "  MIN_CELLS_PER_UNIT:     $MIN_CELLS_PER_UNIT"
echo "  MIN_UNITS_PER_BIN:      $MIN_UNITS_PER_BIN"
echo "  SEED:                   $SEED"
echo "  start:                  $(date -Iseconds)"
echo "=========================================================================="

"$RSCRIPT" "$SCRIPT" \
    --input "$INPUT" \
    --lineage all \
    --pseudotime_col "$PSEUDOTIME_COL" \
    --covariate_col "$COVARIATE_COL" \
    --output_dir "$OUTDIR" \
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
echo "=========================================================================="
echo "Pseudobulk pseudotime run complete: $(date -Iseconds)"
echo "  Run the QC notebook to compare against the cell-level LRT:"
echo "    $SCRIPT_DIR/260526-pseudobulk_concordance.ipynb"
echo "=========================================================================="
