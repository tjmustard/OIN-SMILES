#!/usr/bin/env bash
# v0.4.16 Lane 1: LIVE confirmation of the derived knee, plus the bound-0 wiring gate.
#
# WHY THIS EXISTS EVEN THOUGH THE CURVE IS ALREADY DERIVED. A derived curve is a PREDICTION.
# Reading one as an end-to-end result is the hole v0.4.13 and v0.4.14 both fell into with the
# offline re-score, and the rule that came out of it is that the shipped default is decided by a
# live arm. The derivation says WHICH bound to test; this says whether it holds.
#
# THREE ARMS, and the first one is the instrument check rather than a measurement:
#
#   bound0    OIN_STRING_EXACT_BOUND=0 -- must be BYTE-IDENTICAL to the lever-OFF arm, because
#             bound 0 returns the incumbent the instant it is recorded, which is exactly the
#             conformer the pre-lever code stopped on. 🔴 A NON-ZERO GAIN OR LOSS HERE IS A
#             WIRING BUG, NOT A FINDING, and it invalidates every derived number.
#   target    the chosen bound N over all 365 key_equal -- the accuracy and runtime numbers
#   control   the same bound over 200 byte_exact/INTACT molecules -- the regression arm. A
#             targeted A/B is blind to damage among the 3858 without it.
#
# ⚠ THE BOUND IS INERT WHEN THE LEVER IS OFF, which is what makes this a clean single-variable
# A/B. `incumbent_hit` is only ever set when accept_fn returns ACCEPT_INCUMBENT, and that only
# happens under OIN_ACCEPT_STRING_EXACT -- so exporting the bound into BOTH arms cannot touch the
# OFF arm. Pinned by test_string_exact_bound.py::TestTheBoundCannotFireWithoutAnIncumbent.
#
# ⚠ RUN FROM THE CHECKOUT UNDER TEST. generator_ab_honest.py does sys.path.insert(0, ../src)
# relative to ITSELF, which BEATS PYTHONPATH. Running main's copy silently measures main and
# prints a flawless null -- that cost v0.4.15 six arms.
#
# Usage:  bash tools/run_v0416_confirm.sh <bound_N>
set -uo pipefail

BOUND="${1:?usage: run_v0416_confirm.sh <bound_N>   (the knee from tools/run_v0416_knee.sh)}"
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
V=$MAIN/.venv/bin/python                 # rdkit pinned ==2025.9.3 -- never `uv sync` in a worktree
COHORT=$MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v0.4.5-5k
POP=$MAIN/measurements/v0.4.15
OUT=$MAIN/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.16-confirm

mkdir -p "$OUT"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

dangling=$(find "$COHORT" -xtype l | wc -l)
if [ "$dangling" -ne 0 ]; then
  echo "🔴 REFUSING: $dangling dangling symlinks in $COHORT -- restore the dataset first" >&2
  exit 1
fi

run_arm() {  # name  bound  population_file
  local name=$1 bound=$2 pop=$3
  echo "=== arm $name (bound=$bound) over $(basename "$pop") ==="
  OIN_STRING_EXACT_BOUND="$bound" "$V" "$HERE/tools/generator_ab_honest.py" \
      --cohort-dir "$COHORT" \
      --molecules-file "$pop" \
      --lever OIN_ACCEPT_STRING_EXACT \
      --timeout 300 \
      --out-json "$OUT/ab_bound_${name}.json" > "$OUT/${name}.log" 2>&1 &
}

# The wiring gate runs FIRST and alone. If it fails there is no point paying for the other two.
run_arm bound0 0 "$POP/pop_L2_target_key_equal.txt"
wait
"$V" - "$OUT/ab_bound_bound0.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
g, l = len(d["gains"]), len(d["losses"])
print(f"\n=== BOUND-0 WIRING GATE: gains={g} losses={l} over n={d['n_scored']}")
if g or l:
    print("🔴 FAILED. bound=0 must reproduce the lever-OFF answer byte-for-byte. This is a WIRING")
    print("   BUG, not a finding, and every derived number in this release is void until fixed.")
    raise SystemExit(1)
print("✅ bound=0 is byte-identical to lever-OFF. The zero point of the curve sits where it must.")
PY
[ $? -eq 0 ] || exit 1

run_arm target  "$BOUND" "$POP/pop_L2_target_key_equal.txt"
run_arm control "$BOUND" "$POP/pop_control_byte_exact_intact_200.txt"
wait

echo "=== confirmation arms complete ==="
for a in bound0 target control; do
  "$V" -c "
import json,sys
d=json.load(open(sys.argv[1]))
print(f\"  {sys.argv[2]:<8} n={d['n_scored']:<4} gains={len(d['gains']):<4} losses={len(d['losses'])}\")
" "$OUT/ab_bound_$a.json" "$a" 2>/dev/null || echo "  $a: MISSING"
done
