#!/bin/bash
# Parameterized round-trip sweep runner.
#
# Replaces the one-off run_*_sweep.sh scripts that previous releases dropped into their
# (gitignored) results directories, where they were lost. Same shape as
# results-v0.4.4-regression/run_regression_sweep.sh, but version-controlled and
# parameterized, so a sweep is reproducible from the repo alone.
#
#   tools/run_sweep.sh <cohort-dir> <output-dir> [n_shards] [mol_timeout_s]
#
# Notes that matter for comparability:
#   * <cohort-dir> should be a FROZEN symlink dir from tools/build_sweep_cohort.py, not the
#     raw cat/photo tree. 1,033 basenames exist in BOTH subdirs and the harness keys reports
#     by basename, so the raw tree double-matches them and races their report writes.
#   * mol_timeout is recorded in every report. Do NOT compare pass rates across sweeps with
#     different budgets: v0.4.4's 11 apparent "regressions" were all 300s timeouts against a
#     v0.4.0 arm run under quick mode's 30s cap, i.e. a config artifact, zero wrong answers.
#   * tools/injectivity/missed_success_audit.py hardcodes TIMEOUT_S = 300.0, so a sweep at a
#     different budget will be mis-attributed by it.
#
# Uses the repo's own .venv (rdkit pinned 2025.9.3). Never `uv sync` here.
set -u

COHORT="${1:?usage: run_sweep.sh <cohort-dir> <output-dir> [n_shards] [mol_timeout_s]}"
OUT="${2:?usage: run_sweep.sh <cohort-dir> <output-dir> [n_shards] [mol_timeout_s]}"
N="${3:-6}"
MOL_TIMEOUT="${4:-300}"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="$REPO/.venv/bin/python"

[ -d "$COHORT" ] || { echo "error: cohort dir $COHORT not found" >&2; exit 1; }
mkdir -p "$OUT"
cd "$REPO" || exit 1

NMOL=$(find "$COHORT" -name '*.xyz' | wc -l)
echo "[launch] $(date -Is) HEAD=$(git rev-parse --short HEAD) shards=$N mol_timeout=${MOL_TIMEOUT}s molecules=$NMOL"
echo "[launch] cohort=$COHORT"
echo "[launch] output=$OUT"

# Record the exact configuration alongside the results so the run is self-describing.
{
  echo "{"
  echo "  \"launched_at\": \"$(date -Is)\","
  echo "  \"commit_id\": \"$(git rev-parse --short HEAD)\","
  echo "  \"cohort_dir\": \"$COHORT\","
  echo "  \"molecules\": $NMOL,"
  echo "  \"shards\": $N,"
  echo "  \"mol_timeout_s\": $MOL_TIMEOUT,"
  echo "  \"levers\": {"
  env | grep '^OIN_' | sed 's/\(^[^=]*\)=\(.*\)/    "\1": "\2",/' | sed '$ s/,$//'
  echo "  }"
  echo "}"
} > "$OUT/run_config.json"

pids=()
for i in $(seq 1 "$N"); do
    "$PY" tools/test_dataset_roundtrip.py \
        --dataset-dir "$COHORT" \
        --output-dir "$OUT" \
        --shard "$i:$N" \
        --no-summary \
        --mol-timeout "$MOL_TIMEOUT" \
        > "$OUT/sweep_shard_${i}.log" 2>&1 &
    pids+=($!)
    echo "[launch] shard $i/$N pid=${pids[-1]}"
done

fail=0
for idx in "${!pids[@]}"; do
    wait "${pids[$idx]}" || { echo "[warn] shard $((idx+1)) exited non-zero"; fail=1; }
done
echo "[done] $(date -Is) all shards finished (fail=$fail)"

echo "[merge] rebuilding summary_roundtrip.json from ALL individual_reports/"
"$PY" tools/rebuild_summary.py --output-dir "$OUT"

echo "[report] bucket classification"
PYTHONPATH="$REPO/src" "$PY" tools/roundtrip_bucket_report.py --results-dir "$OUT"

echo "[all-done] $(date -Is)"
