#!/usr/bin/env bash
# ARM 2 golden re-freeze, step 2: run the gate itself. Detached, 6 shards, both goldens.
#
#   run_arm2_refreeze.sh on              the FULL gate, shipped defaults (no OIN_* set), 425 rows
#   run_arm2_refreeze.sh off             OIN_EXACT_DONOR_FOLD=0, ONLY the rows `on` could not
#                                        reproduce (tools/v0417/arm2_refreeze.py diff writes them)
#
# WHY THE FULL GATE AND NOT "THE MOVERS". v0.4.13 and v0.4.14 each re-ran only the rows a sweep
# predicted would move. tools/v0417/arm2_field2_audit.py found rows that were ALREADY stale when
# this lane started -- left behind by exactly that procedure. A predicted list cannot see a row
# nobody predicted, and ARM 2 is only ever run in full at a release, so the gate can sit red
# unseen. One full run is the only instrument that says which rows a re-freeze owes.
#
# WHY `off` AT ALL. A row that differs from the golden is not thereby a row this lever moved.
# Generation is SEEDED (MetalloGenAdapter seed=42), so with the lever off the old row must come
# back byte-for-byte; where it does not, the golden was stale before v0.4.17 and the recorded
# reason has to say so.
#
# The gate exits non-zero on a mismatch. Here that is the expected result, not a failure: the
# rows land in --out either way, and `#DONE n` is what says a shard finished.
set -euo pipefail

ARM="${1:?usage: $0 on|off}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset
OUT="$DATA/results-v0.4.17-exactfold/arm2_refreeze"
SHARDS=6

case "$ARM" in
  on)  LEVER_ENV=(); C047="$DATA/cohort-v047-slow100"; C049="$DATA/cohort-v049-strata" ;;
  off) LEVER_ENV=(-E OIN_EXACT_DONOR_FOLD=0)
       C047="$DATA/cohort-v0.4.17-arm2-control-v047"; C049="$DATA/cohort-v0.4.17-arm2-control-v049" ;;
  *)   echo "usage: $0 on|off" >&2; exit 2 ;;
esac

if [ -n "$(git -C "$HERE" status --porcelain)" ]; then
  echo "REFUSING: $HERE is dirty. The run is attributed to a commit; make it one." >&2; exit 1
fi
if systemctl --user show-environment | grep -q '^OIN_'; then
  echo "REFUSING: the user manager exports an OIN_* variable; 'shipped defaults' would be a lie." >&2
  exit 1
fi
if [ -d "$OUT/$ARM" ] && ls "$OUT/$ARM"/*.tsv >/dev/null 2>&1; then
  echo "REFUSING: $OUT/$ARM already holds results. Move it aside; never append to a run." >&2; exit 1
fi
# The control cohort is a SUBSET and the gate refuses an empty shard, so shard only what exists.
for c in "$C047" "$C049"; do
  [ -d "$c" ] || { echo "REFUSING: cohort $c does not exist" >&2; exit 1; }
  [ "$(find "$c" -xtype l | wc -l)" -eq 0 ] || { echo "REFUSING: dangling symlinks in $c" >&2; exit 1; }
done
mkdir -p "$OUT/$ARM"
git -C "$HERE" log -1 --format='%H %s' > "$OUT/$ARM/COMMIT"

for i in $(seq 1 $SHARDS); do
  cmd=""
  for g in v047:"$C047" v049:"$C049"; do
    tag="${g%%:*}"; cohort="${g#*:}"
    n="$(find "$cohort" -maxdepth 1 -name '*.xyz' | wc -l)"
    # a subset smaller than the shard count leaves some shards empty; those are skipped, loudly
    if [ "$n" -lt "$i" ]; then
      cmd+="echo 'shard $i: $tag has only $n rows, nothing to run' > '$OUT/$ARM/${tag}_$i.log'; "
      cmd+="printf '#DONE 0\n' > '$OUT/$ARM/${tag}_$i.tsv'; "
      continue
    fi
    cmd+="bash '$HERE/tools/gate_v047.sh' arm2 --cohort-dir '$cohort' "
    cmd+="--golden '$HERE/tools/gate_${tag}_arm2_golden.tsv' --shard $i:$SHARDS --timeout 300 "
    cmd+="--out '$OUT/$ARM/${tag}_$i.tsv' > '$OUT/$ARM/${tag}_$i.log' 2>&1 || true; "
  done
  cmd+="date -Is > '$OUT/$ARM/DONE_$i'"
  systemd-run --user --unit="oin-v0417-arm2-$ARM-$i" \
    --description="v0.4.17 ARM 2 re-freeze, arm $ARM, shard $i/$SHARDS" \
    -p OOMPolicy=continue -p MemoryMax=14G \
    "${LEVER_ENV[@]}" -E PATH="$PATH" \
    -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E NUMEXPR_NUM_THREADS=1 \
    /bin/bash -c "$cmd"
done
echo "launched oin-v0417-arm2-$ARM-{1..$SHARDS}. Finished when $OUT/$ARM/DONE_{1..$SHARDS} all exist."
