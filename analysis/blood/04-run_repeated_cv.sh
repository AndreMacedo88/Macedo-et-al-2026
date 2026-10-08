#!/bin/bash
# Launch N stratified-10-fold partitions in parallel (each an independent process with
# CV_WORKERS BiocParallel workers), then aggregate into the final pooled result + plots.
# Idempotent: a partition whose rds already exists is skipped (resumable).
#   usage: run_repeated_cv.sh [N_REPEATS=20] [PAR=8]   (env CV_WORKERS default 4)
set -uo pipefail
: "${M26_DATA_ROOT:?Set M26_DATA_ROOT to the data directory (see README)}"
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="$HERE/04-260702-CV10fold_genelevel_curated_sexadj_ridge.R"
RS="${RSCRIPT:-Rscript}"
PARTDIR="${M26_DATA_ROOT}/3utr/Blood/2025_analyses/CV10fold_genelevel_curated_sexadj_ridge/partitions"
LOGDIR="${M26_RESULTS_ROOT:-$HERE/../../results}/logs/cv10_logs"
mkdir -p "$PARTDIR" "$LOGDIR"
N="${1:-20}"
PAR="${2:-8}"
export CV_WORKERS="${CV_WORKERS:-4}"
BASE_SEED=123

echo "[driver] launching up to $N partitions | PAR=$PAR | CV_WORKERS=$CV_WORKERS"
seq 1 "$N" | xargs -P "$PAR" -I {} bash -c '
  rep="{}"; seed=$(('"$BASE_SEED"' + rep - 1))
  padded=$(printf "%02d" "$rep")
  rds="'"$PARTDIR"'/rep_${padded}_seed_${seed}.rds"
  if [ -f "$rds" ]; then echo "[driver] rep $rep SKIP (exists)"; exit 0; fi
  if '"$RS"' --vanilla "'"$SCRIPT"'" partition "$rep" "$seed" > "'"$LOGDIR"'/rep_${padded}.log" 2>&1; then
    echo "[driver] rep $rep (seed $seed) OK -> $(tail -1 "'"$LOGDIR"'/rep_${padded}.log")"
  else
    echo "[driver] rep $rep (seed $seed) FAILED (see '"$LOGDIR"'/rep_${padded}.log)"
  fi
'
NDONE=$(ls "$PARTDIR"/rep_*.rds 2>/dev/null | wc -l)
echo "[driver] partitions complete: $NDONE / $N present"
echo "[driver] === aggregating ==="
"$RS" --vanilla "$SCRIPT" aggregate "$N" 2>&1
echo "[driver] ALL DONE"
