#!/usr/bin/env bash
# v0.4.18 SELECTION lane -- OIN_ETA_RETARGET through the real sweep harness, over all 1,146 eta-bound
# molecules. ONE arm is run; the OFF arm is the sweep of record's own rows.
#
# WHY NOT TWO ARMS. Shipped defaults ARE results-v0.4.18-sweep, and generation is deterministic at
# this load: the L2 A/B's OFF arm reproduced the v0.4.17 sweep on 1,063/1,063 structures, the
# v0.4.18 sweep reproduced that A/B's ON arm on 1,068/1,068, and this lane's lever-OFF probe
# reproduced the v0.4.18 sweep on 103/103. So ab_off is BUILT, not run: symlinks to that sweep's
# individual reports and structures for the cohort. ⚠ That makes the report's "noise floor" line
# circular here (OFF *is* the record) -- it prints 0 by construction and measures nothing; the
# determinism evidence is the three numbers above.
#
# The arm runs 6 shards, the same six processes on the box as the sweep and both earlier A/Bs.
#
#   tools/v0418/run_selection_ab.sh          # then: AB_OUT=... tools/v0418/post_eta_ab.sh retarget
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
DATA=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset
SWEEP=$DATA/results-v0.4.18-sweep
COHORT=$DATA/cohort-v0.4.18-eta
export AB_OUT=$DATA/results-v0.4.18-eta-selection

[ -f "$SWEEP/DONE" ] || { echo "ABORT: $SWEEP is not a finished sweep"; exit 1; }
off=$AB_OUT/ab_off
if [ ! -f "$off/DONE" ]; then
  rm -rf "$off"; mkdir -p "$off/individual_reports" "$off/structures"
  n=0
  for x in "$COHORT"/*.xyz; do
    m=$(basename "$x" .xyz)
    [ -f "$SWEEP/individual_reports/$m.json" ] || { echo "ABORT: the sweep has no report for $m"; exit 1; }
    ln -s "$SWEEP/individual_reports/$m.json" "$off/individual_reports/$m.json"
    for s in "$SWEEP/structures/$m"[_.]*; do [ -e "$s" ] && ln -s "$s" "$off/structures/$(basename "$s")"; done
    n=$((n + 1))
  done
  echo "results-v0.4.18-sweep @ $(grep -o '"commit_id": "[^"]*"' "$SWEEP/run_config.json") -- the sweep of record's own rows, not a run" > "$off/AB_COMMIT"
  echo "#DONE $n" > "$off/DONE"
  echo "ab_off: $n molecules linked from $SWEEP"
fi
exec "$HERE/tools/v0418/run_eta_ab.sh" retarget
