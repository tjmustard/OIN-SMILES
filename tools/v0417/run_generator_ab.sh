#!/usr/bin/env bash
# v0.4.17 L1b -- OIN_EXACT_DONOR_FOLD through the REAL sweep harness, two arms, one mover cohort.
#
# WHY THE HARNESS AND NOT tools/generator_ab_honest.py. That tool keeps only a byte_exact bit per
# arm: no strings, no generated structure, no SIGKILL budget, and it writes its JSON once at the
# end. This lane owes BOTH numbers -- self-consistent AND verified -- and the verified one is the
# census ruler judging the GENERATED STRUCTURE against the input, so the structures must survive.
# tools/test_dataset_roundtrip.py is what produced the sweep of record (results-v0.4.14-sweep):
# per-molecule killable subprocess, tier ladder, individual reports, structures/*_generated.xyz.
# Same instrument, same flags, same thread caps -- only the lever differs between the arms.
#
# WHY A FRESH OFF ARM instead of reading the sweep of record. OIN3DGenerator(timeout=) is
# advisory, so pool size -- hence accuracy -- depends on machine load. Both arms run AT THE SAME
# TIME, three shards each, so whatever the box does it does to both.
#
#   tools/v0417/run_generator_ab.sh <molecules-file>
#
# A --user unit starts with a minimal PATH; the sweep of record ran from a login shell with
# ~/bin (gxtb) on it, so PATH is passed through. CHECK the `Env:` line of a shard log after
# launch: xtb_available must read True in BOTH arms or the tier ladder is not the sweep's.
#
# Waits on nothing. Completeness = REPORT COUNT per arm == cohort size, never exit status.
set -euo pipefail

MOLS=${1:?molecules file, one name per line}
HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
SRC_COHORT=$DATA/cohort-v0.4.5-5k
COHORT=$DATA/cohort-v0.4.17-exactfold-movers
OUT=$DATA/results-v0.4.17-exactfold
SHARDS=3

[ "$(find "$SRC_COHORT" -xtype l | wc -l)" = 0 ] || { echo "ABORT: dangling cohort symlinks -- restore the dataset"; exit 1; }

mkdir -p "$COHORT"
n=0
while read -r m; do
  [ -z "$m" ] && continue
  case "$m" in \#*) continue ;; esac
  src=$(readlink -f "$SRC_COHORT/$m.xyz")
  [ -f "$src" ] || { echo "ABORT: $m has no input file"; exit 1; }
  ln -sfn "$src" "$COHORT/$m.xyz"
  n=$((n + 1))
done < "$MOLS"
have=$(find "$COHORT" -name '*.xyz' | wc -l)
[ "$have" = "$n" ] || { echo "ABORT: cohort dir has $have files, list has $n -- stale cohort dir?"; exit 1; }
echo "cohort: $n molecules -> $COHORT"

for arm in off on; do
  lever=0; [ "$arm" = on ] && lever=1
  mkdir -p "$OUT/ab_$arm"
  cmd="cd $HERE && for i in \$(seq 1 $SHARDS); do $PY tools/test_dataset_roundtrip.py --dataset-dir $COHORT --output-dir $OUT/ab_$arm --shard \$i:$SHARDS --mol-timeout 300 --no-summary > $OUT/ab_$arm/shard\$i.log 2>&1 & done; wait; echo \"#DONE \$(ls $OUT/ab_$arm/individual_reports | wc -l)\" > $OUT/ab_$arm/DONE"
  systemd-run --user --unit="oin-v0417-ab-$arm" \
    --description="v0.4.17 L1b generator A/B, arm $arm (OIN_EXACT_DONOR_FOLD=$lever)" \
    -p OOMPolicy=continue -p MemoryMax=12G \
    -E OIN_EXACT_DONOR_FOLD=$lever -E PYTHONPATH="$HERE/src" -E PATH="$PATH" \
    -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E NUMEXPR_NUM_THREADS=1 \
    /bin/bash -c "$cmd"
done
echo "launched: oin-v0417-ab-off, oin-v0417-ab-on  ($SHARDS shards each). Done when both $OUT/ab_*/DONE exist."
