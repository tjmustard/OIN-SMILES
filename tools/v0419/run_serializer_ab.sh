#!/usr/bin/env bash
# v0.4.19 serializer lane -- the two string-side levers through the REAL sweep harness, over the
# CHANGED-STRING SET: every molecule whose input string moves under
#
#   OIN_H_FAITHFUL=1  OIN_RC1_PROPAGATE=1
#
# (the set comes from the full-cohort e_selfconsistency run, base column vs the release sweep's
# smiles_1 -- tools/v0419/serializer_ab_report.py --changed-set writes it). A molecule whose string
# does NOT move is unchanged by construction on the generator side (generation is byte-deterministic
# at this load: v0.4.18's three reproductions, 1,063 + 1,068 + 103 of as many) and is re-scored
# OFFLINE instead (tools/honest_rescore.py under the levers), so it is not in this cohort.
#
# ONE arm is run. The OFF arm IS the sweep of record (shipped defaults, same commit lineage, same
# strings): ab_off is BUILT from symlinks to results-v0.4.18-release-sweep, not re-run -- which makes
# a "noise floor" line circular here; the determinism evidence is the three numbers above.
#
#   tools/v0419/run_serializer_ab.sh <changed_set.txt> [fix2|cap|fix3]   # then: post_serializer_ab.sh [arm]
#     fix2  OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1  (default; ab/)   cap  OIN_CAP_IGNORES_METAL=1 (ab_cap/)   fix3  all three (ab_fix3/)
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
SRC_COHORT=$DATA/cohort-v0.4.5-5k
SWEEP=$DATA/results-v0.4.18-release-sweep
LIST=${1:?usage: run_serializer_ab.sh <changed_set.txt> [fix2|cap|fix3]}
ARM=${2:-fix2}
case "$ARM" in
  fix2) LEVERS="-E OIN_H_FAITHFUL=1 -E OIN_RC1_PROPAGATE=1"; SUB=ab;;
  cap)  LEVERS="-E OIN_CAP_IGNORES_METAL=1"; SUB=ab_cap;;
  fix3) LEVERS="-E OIN_H_FAITHFUL=1 -E OIN_RC1_PROPAGATE=1 -E OIN_CAP_IGNORES_METAL=1"; SUB=ab_fix3;;
  *) echo "unknown arm $ARM"; exit 1;;
esac
OUT=${AB_OUT:-$DATA/results-v0.4.19-serializer/$SUB}
COHORT=$DATA/cohort-v0.4.19-changed-$ARM
SHARDS=${SHARDS:-6}

[ -z "$(git -C "$HERE" status --porcelain)" ] || { echo "REFUSING: $HERE is dirty. The run is attributed to a commit; make it one."; exit 1; }
if systemctl --user show-environment | grep -q '^OIN_'; then echo "REFUSING: the user manager exports an OIN_* variable"; exit 1; fi
[ "$(find "$SRC_COHORT" -xtype l | wc -l)" = 0 ] || { echo "ABORT: dangling cohort symlinks -- restore the dataset"; exit 1; }
[ -f "$SWEEP/DONE" ] || { echo "ABORT: $SWEEP is not a finished sweep"; exit 1; }
[ ! -d "$OUT/ab_on/individual_reports" ] || { echo "REFUSING: $OUT/ab_on already holds reports. Move it aside."; exit 1; }

mkdir -p "$COHORT" "$OUT"
find "$COHORT" -maxdepth 1 -type l -delete     # a stale member would be run as if it belonged
n=0
while read -r m; do
  [ -n "$m" ] || continue; case "$m" in \#*) continue;; esac
  src=$(readlink -f "$SRC_COHORT/$m.xyz")
  [ -f "$src" ] || { echo "ABORT: $m has no input file"; exit 1; }
  ln -sfn "$src" "$COHORT/$m.xyz"; n=$((n + 1))
done < "$LIST"
have=$(find "$COHORT" -name '*.xyz' | wc -l)
[ "$have" = "$n" ] && [ "$n" -gt 0 ] || { echo "ABORT: cohort dir has $have files, list says $n"; exit 1; }
echo "cohort: $n changed-string molecules -> $COHORT"

off=$OUT/ab_off
if [ ! -f "$off/DONE" ]; then
  rm -rf "$off"; mkdir -p "$off/individual_reports" "$off/structures"
  for x in "$COHORT"/*.xyz; do
    m=$(basename "$x" .xyz)
    [ -f "$SWEEP/individual_reports/$m.json" ] || { echo "ABORT: the sweep has no report for $m"; exit 1; }
    ln -s "$SWEEP/individual_reports/$m.json" "$off/individual_reports/$m.json"
    for s in "$SWEEP/structures/$m"[_.]*; do [ -e "$s" ] && ln -s "$s" "$off/structures/$(basename "$s")"; done
  done
  echo "results-v0.4.18-release-sweep @ $(grep -o '"commit_id": "[^"]*"' "$SWEEP/run_config.json") -- the sweep of record's own rows, not a run" > "$off/AB_COMMIT"
  echo "#DONE $n" > "$off/DONE"
  echo "ab_off: $n molecules linked from $SWEEP"
fi

mkdir -p "$OUT/ab_on"
git -C "$HERE" log -1 --format='%H %s' > "$OUT/ab_on/AB_COMMIT"
cmd="cd $HERE && for i in \$(seq 1 $SHARDS); do $PY tools/test_dataset_roundtrip.py --dataset-dir $COHORT --output-dir $OUT/ab_on --shard \$i:$SHARDS --mol-timeout 300 --no-summary > $OUT/ab_on/shard\$i.log 2>&1 & done; wait; echo \"#DONE \$(ls $OUT/ab_on/individual_reports | wc -l)\" > $OUT/ab_on/DONE"
systemd-run --user --unit="oin-v0419-serab-$ARM" \
  --description="v0.4.19 A/B, ON arm $ARM ($LEVERS) over its changed-string set" \
  -p OOMPolicy=continue -p MemoryMax=12G \
  $LEVERS \
  -E PYTHONPATH="$HERE/src" -E PATH="$PATH" \
  -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E NUMEXPR_NUM_THREADS=1 \
  /bin/bash -c "$cmd"
echo "launched: ab_on ($SHARDS shards). Done when $OUT/ab_on/DONE exists."
