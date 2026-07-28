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

# Locate the interpreter. A git WORKTREE has no .venv of its own -- this project's rule is
# to use the MAIN checkout's pinned venv and never `uv sync` in a worktree -- and sweeping
# from a worktree is the normal case here, so $REPO/.venv is usually absent.
#
# 🔴 v0.4.13: THE SIBLING GLOB IS GONE. The previous fallback was
#     for cand in "$(dirname "$REPO")"/*/.venv/bin/python; do ... break; done
# which matches EVERY project sharing the parent directory and takes the alphabetically
# first one. gate_v047.sh was hardened against exactly this in v0.4.9 and this script was
# not, so the trap stayed live in the one place where it is most expensive: a 5000-molecule
# sweep is ~55 CPU-h, and an interpreter with the wrong rdkit does not fail -- it produces a
# full, plausible, INCOMPARABLE table. Measured from the v0.4.13 worktree, the glob's first
# match was `EtaCatalysis/.venv` (no rdkit at all) and its second was `EtaTMCSMILES/.venv`
# (rdkit 2025.09.2 against this project's pinned 2025.9.3) -- the loud failure and the silent
# one, in that order.
#
# A worktree resolves to its own main checkout deterministically via `--git-common-dir`, so
# there is nothing to guess. Resolution order: explicit override, this repo's own venv, the
# main checkout's venv. Then the rdkit version is REQUIRED to match the pyproject pin: a
# sweep on a different rdkit cannot be quoted beside any earlier table.
_rdkit_version() { "$1" -c 'import rdkit; print(rdkit.__version__)' 2>/dev/null; }

PIN="$(sed -n 's/.*"rdkit==\([^"]*\)".*/\1/p' "$REPO/pyproject.toml" | head -1)"
MAIN_CHECKOUT="$(dirname "$(git -C "$REPO" rev-parse --git-common-dir 2>/dev/null || echo /nonexistent/.git)")"
PY=""
_tried=""
for cand in "${OIN_SWEEP_PYTHON:-}" "$REPO/.venv/bin/python" "$MAIN_CHECKOUT/.venv/bin/python"; do
    [ -z "$cand" ] || [ ! -x "$cand" ] && { _tried="$_tried ${cand:-<unset>}"; continue; }
    _tried="$_tried $cand"
    [ -n "$(_rdkit_version "$cand")" ] && { PY="$cand"; break; }
done
if [ -z "$PY" ]; then
    echo "error: no python with an importable rdkit found. Set OIN_SWEEP_PYTHON to the MAIN" >&2
    echo "       checkout's .venv/bin/python (rdkit is pinned ==$PIN there; never 'uv sync'" >&2
    echo "       in a worktree). Tried:$_tried" >&2
    exit 1
fi
RDKIT_V="$(_rdkit_version "$PY")"
# Normalize 2025.09.2 vs 2025.9.2 before comparing -- zero-padding differs between the
# pyproject pin and rdkit.__version__ and a spurious warning trains people to ignore it.
_norm() { echo "$1" | awk -F. '{for(i=1;i<=NF;i++){printf "%s%d", (i>1?".":""), $i}; print ""}'; }
if [ -n "$PIN" ] && [ "$(_norm "$RDKIT_V")" != "$(_norm "$PIN")" ]; then
    echo "error: rdkit $RDKIT_V != pinned $PIN. A sweep on a different rdkit produces a table" >&2
    echo "       that cannot be compared to any earlier one. Refusing." >&2
    exit 1
fi
echo "[launch] interpreter=$PY (rdkit $RDKIT_V, pin $PIN)"

[ -d "$COHORT" ] || { echo "error: cohort dir $COHORT not found" >&2; exit 1; }
mkdir -p "$OUT"
cd "$REPO" || exit 1

# The harness appends its OWN ../src to sys.path, so running a copy from a different
# checkout than the code under test silently mixes two trees. Pin it explicitly.
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"
echo "[launch] PYTHONPATH=$PYTHONPATH"

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
