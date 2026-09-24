#!/usr/bin/env bash
# v0.4.19 serializer lane -- after run_serializer_ab.sh's ON arm has written DONE: bucket both arms
# honestly, judge each arm's generated structures with the neutral ruler, parse each arm's OWN
# strings back (the ON arm's smiles_1 is the lever's -- the string clause of VERIFIED must be
# recomputed per arm, unlike v0.4.18's generator-side A/B), then the report.
#
#   tools/v0419/post_serializer_ab.sh
#
# Refuses an incomplete arm: completeness is REPORT COUNT == cohort size.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
COHORT=$DATA/cohort-v0.4.19-changed
LANE=$DATA/results-v0.4.19-serializer
OUT=${AB_OUT:-$LANE/ab}
export PYTHONPATH=$HERE/src

n=$(find "$COHORT" -name '*.xyz' | wc -l)
[ "$n" -gt 0 ] || { echo "ABORT: $COHORT is empty"; exit 1; }
for arm in off on; do
  d=$OUT/ab_$arm
  [ -f "$d/DONE" ] || { echo "ABORT: $d/DONE missing -- arm $arm still running (systemctl --user status oin-v0419-serab-$arm)"; exit 1; }
  r=$(find "$d/individual_reports" -name '*.json' | wc -l)
  [ "$r" = "$n" ] || { echo "ABORT: arm $arm has $r reports for a cohort of $n -- a shard died; re-run it with --continue UNDER THE SAME LEVERS"; exit 1; }
done

cd "$HERE"
for arm in off on; do
  d=$OUT/ab_$arm
  [ -f "$d/bucket_report_honest.json" ] || $PY tools/roundtrip_bucket_report.py --results-dir "$d" --score honest > "$d/bucket_report_honest.log" 2>&1
  [ -f "$d/g_verdict.jsonl" ] || $PY tools/census/g_vs_input.py --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 8 > "$d/g_verdict.log" 2>&1
  [ -f "$d/parseback.jsonl" ] || $PY tools/census/string_sufficiency.py parseback --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 8 > "$d/parseback.log" 2>&1
  [ -f "$d/parseback_gen.jsonl" ] || $PY tools/census/string_sufficiency.py parseback --side gen --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 8 > "$d/parseback_gen.log" 2>&1
done
$PY tools/v0419/serializer_ab_report.py all --ab "$OUT" --out "$LANE/serializer_ab_report.json" | tee "$LANE/serializer_ab_report.txt"
