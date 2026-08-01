#!/usr/bin/env bash
# v0.4.16 Lane 1: ONE instrumented unbounded run over the 365 key_equal molecules, from which
# BOTH the recovered-vs-bound and runtime-vs-bound curves are derived arithmetically.
#
# WHY ONE RUN AND NOT A PARAMETER SWEEP. The v0.4.16 charter budgets "~1-2 h per point" varying
# OIN_STRING_EXACT_BOUND. That is unnecessary, and the reason is structural rather than clever:
# `incumbent_hit` is the FIRST ACCEPT_INCUMBENT conformer and `generate_3d_structures` returns it
# as the sole pool member regardless of how much longer the pool fills. So bounding at N changes
# the answer ONLY for a molecule whose string-exact hit lies beyond N -- and the probe records
# each molecule's `min_bound` (the smallest bound that still recovers it) plus an elapsed stamp at
# every post-incumbent evaluation. Both curves are then arithmetic over one run.
#
# ⚠ DERIVED IS NOT MEASURED. The chosen point still gets a live confirmation arm
# (tools/generator_ab_honest.py). Reading a derived curve as an end-to-end result is precisely the
# hole v0.4.13/v0.4.14 fell into with the offline re-score.
#
# ⚠ RUN THIS FROM THE CHECKOUT UNDER TEST. selection_pool_probe.py does sys.path.insert(0, ../src)
# relative to ITSELF, which BEATS PYTHONPATH. Running main's copy silently measures main's code
# and prints a flawless null -- that cost v0.4.15 six arms and an hour of belief.
#
# ⚠ LOAD BIASES ACCURACY, NOT JUST TIMING. OIN3DGenerator(timeout=) is advisory and the embed loop
# checks its deadline BETWEEN attempts, so CPU starvation shrinks the pool -- which here would
# understate the recovered CEILING, not merely slow things down. BLAS is capped to 1 per process
# and concurrency held to 3, matching the conditions the v0.4.15 arm was measured under so the two
# are comparable.
#
# Usage:  bash tools/run_v0416_knee.sh [n_shards]
set -uo pipefail

MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # the checkout under test
V=$MAIN/.venv/bin/python                 # rdkit pinned ==2025.9.3 -- never `uv sync` in a worktree
COHORT=$MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k
POP=$MAIN/measurements/v0.4.15/pop_L2_target_key_equal.txt
OUT=$MAIN/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.16-knee
NSHARD="${1:-3}"

# Never write into a worktree: those files die on `git worktree remove`.
mkdir -p "$OUT"

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# 🔴 The dataset can vanish: a Kulik_TMC_Dataset checkout and back deletes 26,232 files while
# `git status` stays clean, leaving every cohort symlink dangling and every run silently empty.
dangling=$(find "$COHORT" -xtype l | wc -l)
if [ "$dangling" -ne 0 ]; then
  echo "🔴 REFUSING: $dangling dangling symlinks in $COHORT -- restore the dataset first" >&2
  exit 1
fi
[ -s "$POP" ] || { echo "🔴 REFUSING: population file $POP is empty or missing" >&2; exit 1; }

# Split the population into NSHARD interleaved chunks. Interleaved, not blocked: the population is
# alphabetical and cost is heavily molecule-dependent, so contiguous blocks finish at wildly
# different times while round-robin balances them.
# ⚠ `#` comment lines carry the population's provenance and must not be read as molecule names.
grep -v '^[[:space:]]*#' "$POP" | grep -v '^[[:space:]]*$' > "$OUT/pop_all.txt"
total=$(wc -l < "$OUT/pop_all.txt")
for s in $(seq 1 "$NSHARD"); do
  awk -v s="$s" -v n="$NSHARD" 'NR % n == (s % n)' "$OUT/pop_all.txt" > "$OUT/shard$s.txt"
done
echo "population $total -> $NSHARD shards: $(for s in $(seq 1 "$NSHARD"); do wc -l < "$OUT/shard$s.txt" | tr -d ' '; echo -n ' '; done)"

for s in $(seq 1 "$NSHARD"); do
  "$V" "$HERE/tools/selection_pool_probe.py" \
      --molecules-file "$OUT/shard$s.txt" \
      --cohort-dir "$COHORT" \
      --lever OIN_ACCEPT_STRING_EXACT \
      --limit 100000 \
      --timeout 300 \
      --out-json "$OUT/knee_shard$s.json" > "$OUT/shard$s.log" 2>&1 &
done
wait

echo "=== shards complete ==="
for s in $(seq 1 "$NSHARD"); do
  n=$("$V" -c "import json,sys; print(len(json.load(open(sys.argv[1]))['rows']))" "$OUT/knee_shard$s.json" 2>/dev/null || echo "MISSING")
  echo "  shard$s: $n rows"
done
# Completeness is a ROW COUNT, never an exit status -- a shard that dies mid-run still exits 0
# through the pipe and leaves a plausible, short JSON.
echo "expected total: $total"
