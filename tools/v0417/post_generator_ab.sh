#!/usr/bin/env bash
# v0.4.17 L1b -- after BOTH arms of run_generator_ab.sh have written their DONE file:
# bucket each arm honestly, judge each arm's generated structures with the neutral ruler, report.
#
#   tools/v0417/post_generator_ab.sh
#
# Refuses to run on an incomplete arm: completeness is REPORT COUNT == cohort size, because a
# shard that died early exits 0 like one that finished.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
COHORT=$DATA/cohort-v0.4.17-exactfold-movers
OUT=$DATA/results-v0.4.17-exactfold
export PYTHONPATH=$HERE/src

n=$(find "$COHORT" -name '*.xyz' | wc -l)
for arm in off on; do
  d=$OUT/ab_$arm
  [ -f "$d/DONE" ] || { echo "ABORT: $d/DONE missing -- arm $arm still running (systemctl --user status oin-v0417-ab-$arm)"; exit 1; }
  r=$(find "$d/individual_reports" -name '*.json' | wc -l)
  [ "$r" = "$n" ] || { echo "ABORT: arm $arm has $r reports for a cohort of $n -- a shard died; re-run it with --continue"; exit 1; }
done

cd "$HERE"
for arm in off on; do
  d=$OUT/ab_$arm
  $PY tools/roundtrip_bucket_report.py --results-dir "$d" --score honest > "$d/bucket_report_honest.log" 2>&1
  $PY tools/census/g_vs_input.py --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 10 > "$d/g_verdict.log" 2>&1
done
$PY tools/v0417/generator_ab_report.py | tee "$OUT/generator_ab_report.txt"
