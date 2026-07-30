#!/usr/bin/env bash
# v0.4.15 three-arm A/B: L1-only, L2-only, both -- real generation, honest scoring.
#
# WHY THREE ARMS. The owner accepted two headline movers in one release knowing a reader cannot
# attribute a move to either lane without reading both apart. A combined-only number is
# unattributable and must not be quoted as either lane's.
#
# WHY NOT THREE FULL SWEEPS. Three 5k sweeps is ~24 h wall / ~116 CPU-h; owner chose targeted arms
# plus ONE full sweep of the chosen default (2026-07-29). The per-arm attribution therefore carries
# a sampling caveat, which the release note must state explicitly.
#
# ⚠ generator_ab_honest.py re-perceives the WRITTEN XYZ, never res.mol. Do NOT substitute
# gate_v047.sh arm2 -- that is a byte-identity gate scoring with the generator's own bond graph,
# the circular predicate OIN_INDEP_SCORE replaced at a 9.6% false-positive rate.
#
# ⚠ Load biases ACCURACY, not just timing: OIN3DGenerator(timeout=) is advisory and the embed loop
# checks its deadline BETWEEN attempts, so CPU starvation shrinks the pool and UNDERSTATES
# byte_exact. BLAS is capped to 1 per process below and concurrency is held to 3.
#
# Usage:  bash tools/run_v0415_arms.sh <lane1|lane2>
set -uo pipefail

MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
V=$MAIN/.venv/bin/python                 # rdkit pinned ==2025.9.3 -- never `uv sync` in a worktree
COHORT=$MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k
POP=$MAIN/measurements/v0.4.15
OUT=$MAIN/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.15-arms

# Never a worktree: worktree files die on `git worktree remove`.
mkdir -p "$OUT"

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# 🔴 The dataset can vanish: a Kulik_TMC_Dataset checkout and back deletes 26,232 files while
# `git status` stays clean, leaving every cohort symlink dangling.
dangling=$(find "$COHORT" -xtype l | wc -l)
if [ "$dangling" -ne 0 ]; then
  echo "🔴 REFUSING: $dangling dangling symlinks in $COHORT -- restore the dataset first" >&2
  exit 1
fi

WHICH=${1:?usage: run_v0415_arms.sh <lane1|lane2|both>}
case "$WHICH" in
  lane1) SRC=$MAIN/../oin-v0415-attach;     LEVER=OIN_ATTACH_RETURN
         POPS="pop_L1_target_site_lost pop_L1_control_byte_exact_detached
               pop_control_byte_exact_intact_200" ;;
  lane2) SRC=$MAIN/../oin-v0415-enantiomer; LEVER=OIN_ACCEPT_STRING_EXACT
         POPS="pop_L2_mirror_match pop_L2_target_key_equal
               pop_control_byte_exact_intact_200" ;;
  both)  # The combined arm, decomposed the way the lanes land: Lane 1 HELD ON, Lane 2 varying.
         # So this reads as Lane 2's MARGINAL effect on top of Lane 1, which is the number a
         # reader needs to attribute the combined move -- not an unattributable joint total.
         # `_run_one` only sets `--lever` and inherits the rest of the environment, so exporting
         # the holding lever here keeps it on in BOTH sides of the A/B.
         SRC=$MAIN/../oin-v0415-both;       LEVER=OIN_ACCEPT_STRING_EXACT
         export OIN_ATTACH_RETURN=1
         POPS="pop_L2_mirror_match pop_L1_target_site_lost
               pop_control_byte_exact_intact_200" ;;
  *) echo "unknown arm $WHICH" >&2; exit 2 ;;
esac

if [ ! -d "$SRC/src" ]; then
  echo "🔴 REFUSING: $SRC/src does not exist -- create the worktree first" >&2
  exit 1
fi

# 🔴 EACH ARM MUST RUN THE TOOL OUT OF ITS OWN CHECKOUT, and PYTHONPATH IS NOT ENOUGH.
#
# `generator_ab_honest.py` line 78 does
#     sys.path.insert(0, <dirname(__file__)>/../src)
# which puts ITS OWN checkout's src at position 0 and therefore OVERRIDES PYTHONPATH. The first
# version of this script set PYTHONPATH=$SRC/src but invoked $MAIN/tools/generator_ab_honest.py,
# so all six arms imported MAIN's oinsmiles -- where neither v0.4.15 lever exists. Both sides of
# every A/B ran identical code and every arm returned a perfectly clean
# "0 gains, 0 losses, output moved 0". A broken A/B prints exactly the null a real one would.
#
# Invoking $SRC/tools/... turns that self-locating insert from a trap into the guarantee: the
# tool, the code it imports and the lever under test all come from one tree. PYTHONPATH is kept
# as a belt-and-braces second signal, not as the mechanism.
TOOL=$SRC/tools/generator_ab_honest.py
export PYTHONPATH=$SRC/src

# Prove it before spending hours: the resolved package must live under $SRC, not $MAIN.
resolved=$("$V" -c "import sys,os; sys.path.insert(0, os.path.join('$SRC','src')); import oinsmiles; print(oinsmiles.__file__)")
case "$resolved" in
  "$SRC"/src/oinsmiles/*) echo "  code under test: $resolved" ;;
  *) echo "🔴 REFUSING: oinsmiles resolves to $resolved, not $SRC/src -- the arm would measure the wrong tree" >&2
     exit 1 ;;
esac

echo "=== $WHICH: lever $LEVER, code $SRC ==="
for p in $POPS; do
  n=$(grep -cvE '^\s*(#|$)' "$POP/$p.txt")
  echo "--- $p (n=$n) -> $OUT/${WHICH}_${p}.json"
  nohup "$V" "$TOOL" \
      --cohort-dir "$COHORT" \
      --molecules-file "$POP/$p.txt" \
      --lever "$LEVER" \
      --timeout 300 \
      --out-json "$OUT/${WHICH}_${p}.json" \
      > "$OUT/${WHICH}_${p}.log" 2>&1 &
  echo "    pid $!"
done
wait
echo "=== $WHICH DONE ===" | tee -a "$OUT/PROGRESS"
echo "${WHICH}_ARMS_DONE" >> "$OUT/PROGRESS"
