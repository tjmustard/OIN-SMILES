#!/usr/bin/env bash
# v0.4.18 L2 -- the eta construction levers through the REAL sweep harness, two arms, one cohort.
#
#   off   OIN_ETA_COVALENT_TARGET=0  OIN_VDW_EXEMPT_BINDING=0  OIN_ETA_TARGET_UNSCALED=0   (written, never unset)
#   on    OIN_ETA_COVALENT_TARGET=1  OIN_VDW_EXEMPT_BINDING=1  OIN_ETA_TARGET_UNSCALED=0
#   on3   OIN_ETA_COVALENT_TARGET=1  OIN_VDW_EXEMPT_BINDING=1  OIN_ETA_TARGET_UNSCALED=1
#
# THE COHORT IS EVERY ETA-BOUND MOLECULE (1,146 of 5,000), not "the detached ones". Both levers are
# read only inside a complex that HAS a multi-atom binding group -- the clash exemption covers every
# binding atom of such a complex, and nothing in any other -- so a non-eta molecule is byte-identical
# by construction (checked: tools/v0418/eta_path_probe.py on non-eta molecules of every dead-zone
# metal, 25 of 25 built structures identical). The eta molecules that PASS today are in the cohort
# because they are where a loss would come from: the covalent target moves EVERY eta ring.
#
# Same instrument as the sweep of record (tools/test_dataset_roundtrip.py, --mol-timeout 300,
# BLAS=1), and BOTH ARMS AT ONCE, three shards each: OIN3DGenerator(timeout=) is advisory, so
# accuracy depends on load, and whatever the box does it must do to both.
#
#   tools/v0418/run_eta_ab.sh            # arms off + on, 3 shards each, at once
#   tools/v0418/run_eta_ab.sh on3        # ONE more arm, 6 shards -- the same six processes on the box.
#                                        # Licensed by the first run: its OFF arm reproduced the sweep
#                                        # of record byte-for-byte (1,063/1,063), so OFF need not be re-run.
#
# Waits on nothing. Completeness = REPORT COUNT per arm == cohort size, never exit status.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
SRC_COHORT=$DATA/cohort-v0.4.5-5k
COHORT=$DATA/cohort-v0.4.18-eta
OUT=$DATA/results-v0.4.18-eta-detached
TABLE=$DATA/results-v0.4.17-reattribution/attribution_table.tsv
ARMS=("$@"); [ ${#ARMS[@]} -gt 0 ] || ARMS=(off on)
SHARDS=$((6 / ${#ARMS[@]}))
levers() {  # arm -> "COVALENT EXEMPT UNSCALED"
  case "$1" in off) echo "0 0 0";; on) echo "1 1 0";; on3) echo "1 1 1";; *) echo "unknown arm $1" >&2; exit 1;; esac
}
for arm in "${ARMS[@]}"; do levers "$arm" > /dev/null; done

[ -z "$(git -C "$HERE" status --porcelain)" ] || { echo "REFUSING: $HERE is dirty. The run is attributed to a commit; make it one."; exit 1; }
if systemctl --user show-environment | grep -q '^OIN_'; then echo "REFUSING: the user manager exports an OIN_* variable"; exit 1; fi
[ "$(find "$SRC_COHORT" -xtype l | wc -l)" = 0 ] || { echo "ABORT: dangling cohort symlinks -- restore the dataset"; exit 1; }
for arm in "${ARMS[@]}"; do
  [ ! -d "$OUT/ab_$arm/individual_reports" ] || { echo "REFUSING: $OUT/ab_$arm already holds reports. Move it aside."; exit 1; }
done

mkdir -p "$COHORT"
find "$COHORT" -maxdepth 1 -type l -delete     # a stale member would be run as if it belonged
n=0
while read -r m; do
  src=$(readlink -f "$SRC_COHORT/$m.xyz")
  [ -f "$src" ] || { echo "ABORT: $m has no input file"; exit 1; }
  ln -sfn "$src" "$COHORT/$m.xyz"; n=$((n + 1))
done < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $c["eta"]==1 {print $c["molecule"]}' "$TABLE")
have=$(find "$COHORT" -name '*.xyz' | wc -l)
[ "$have" = "$n" ] && [ "$n" -gt 0 ] || { echo "ABORT: cohort dir has $have files, table says $n"; exit 1; }
echo "cohort: $n eta-bound molecules -> $COHORT"
for arm in "${ARMS[@]}"; do
  read -r cov exempt unscaled <<< "$(levers "$arm")"
  mkdir -p "$OUT/ab_$arm"
  git -C "$HERE" log -1 --format='%H %s' > "$OUT/ab_$arm/AB_COMMIT"
  cmd="cd $HERE && for i in \$(seq 1 $SHARDS); do $PY tools/test_dataset_roundtrip.py --dataset-dir $COHORT --output-dir $OUT/ab_$arm --shard \$i:$SHARDS --mol-timeout 300 --no-summary > $OUT/ab_$arm/shard\$i.log 2>&1 & done; wait; echo \"#DONE \$(ls $OUT/ab_$arm/individual_reports | wc -l)\" > $OUT/ab_$arm/DONE"
  systemd-run --user --unit="oin-v0418-etaab-$arm" \
    --description="v0.4.18 L2 eta A/B, arm $arm (COVALENT_TARGET=$cov EXEMPT_BINDING=$exempt TARGET_UNSCALED=$unscaled)" \
    -p OOMPolicy=continue -p MemoryMax=12G \
    -E OIN_ETA_COVALENT_TARGET=$cov -E OIN_VDW_EXEMPT_BINDING=$exempt -E OIN_ETA_TARGET_UNSCALED=$unscaled \
    -E PYTHONPATH="$HERE/src" -E PATH="$PATH" \
    -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E NUMEXPR_NUM_THREADS=1 \
    /bin/bash -c "$cmd"
done
echo "launched: ${ARMS[*]} ($SHARDS shards each). Done when $OUT/ab_<arm>/DONE exists for each."
