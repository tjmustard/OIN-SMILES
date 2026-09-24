#!/usr/bin/env bash
# v0.4.18 L2 -- after BOTH arms of run_eta_ab.sh have written their DONE file: bucket each arm
# honestly, judge each arm's generated structures with the neutral ruler, report both numbers.
#
#   tools/v0418/post_eta_ab.sh [on|on3]     # the ON arm to read against ab_off (default on)
#
# Refuses an incomplete arm: completeness is REPORT COUNT == cohort size, because a shard that
# died early exits 0 like one that finished.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
COHORT=$DATA/cohort-v0.4.18-eta
OUT=$DATA/results-v0.4.18-eta-detached
ON=${1:-on}
export PYTHONPATH=$HERE/src

n=$(find "$COHORT" -name '*.xyz' | wc -l)
[ "$n" -gt 0 ] || { echo "ABORT: $COHORT is empty"; exit 1; }
for arm in off "$ON"; do
  d=$OUT/ab_$arm
  [ -f "$d/DONE" ] || { echo "ABORT: $d/DONE missing -- arm $arm still running (systemctl --user status oin-v0418-etaab-$arm)"; exit 1; }
  r=$(find "$d/individual_reports" -name '*.json' | wc -l)
  [ "$r" = "$n" ] || { echo "ABORT: arm $arm has $r reports for a cohort of $n -- a shard died; re-run it with --continue UNDER THE SAME LEVERS"; exit 1; }
done

cd "$HERE"
for arm in off "$ON"; do
  d=$OUT/ab_$arm
  [ -f "$d/g_verdict.jsonl" ] && [ -f "$d/bucket_report_honest.json" ] && continue   # already read
  $PY tools/roundtrip_bucket_report.py --results-dir "$d" --score honest > "$d/bucket_report_honest.log" 2>&1
  $PY tools/census/g_vs_input.py --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 10 > "$d/g_verdict.log" 2>&1
done
sfx=""; [ "$ON" = on ] || sfx="_$ON"
$PY tools/v0418/eta_ab_report.py --on-arm "$ON" --out "$OUT/eta_ab_report$sfx.json" | tee "$OUT/eta_ab_report$sfx.txt"
