#!/usr/bin/env bash
# v0.4.19 serializer lane -- after run_serializer_ab.sh's ON arm has written DONE: bucket both arms
# honestly, judge each arm's generated structures with the neutral ruler, parse each arm's OWN
# strings back (the ON arm's smiles_1 is the lever's -- the string clause of VERIFIED must be
# recomputed per arm, unlike v0.4.18's generator-side A/B), then the report.
#
#   tools/v0419/post_serializer_ab.sh [fix2|cap|fix3|e2b|e2c|e2f|res]     (e2* and res arms live in results-v0.4.19-e2/)
#
# Refuses an incomplete arm: completeness is REPORT COUNT == cohort size.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
LANE=$DATA/results-v0.4.19-serializer
ARM=${1:-fix2}
case "$ARM" in
  fix2) ON_LEVERS="OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1"; SUB=ab; ESC=e_selfconsistency_fix2.jsonl; RESCORE=rescore_fix2;;
  cap)  ON_LEVERS="OIN_CAP_IGNORES_METAL=1"; SUB=ab_cap; ESC=e_selfconsistency_cap.jsonl; RESCORE=rescore_cap;;
  fix3) ON_LEVERS="OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1 OIN_CAP_IGNORES_METAL=1"; SUB=ab_fix3; ESC=e_selfconsistency_fix3.jsonl; RESCORE=rescore_fix3;;
  e2b)  ON_LEVERS="OIN_N_VALENCE_2=1 OIN_CANONICAL_RESONANCE=1"; SUB=ab_e2b; ESC=e_selfconsistency_e2b.jsonl; RESCORE=rescore_e2b;;
  e2c)  ON_LEVERS="OIN_N_VALENCE_2=1 OIN_CANONICAL_RESONANCE=1 OIN_CANONICAL_CHARGES=1"; SUB=ab_e2c; ESC=e_selfconsistency_e2c.jsonl; RESCORE=rescore_e2c;;
  e2f)  ON_LEVERS="OIN_N_VALENCE_2=1 OIN_CANONICAL_RESONANCE=1"; SUB=ab_e2f; ESC=e_selfconsistency_e2f.jsonl; RESCORE=rescore_e2f;;
  res)  ON_LEVERS="OIN_CANONICAL_RESONANCE=1"; SUB=ab_res; ESC=e_selfconsistency_res.jsonl; RESCORE=rescore_res;;
  *) echo "unknown arm $ARM"; exit 1;;
esac
case "$ARM" in e2*|res) LANE=$DATA/results-v0.4.19-e2;; esac
COHORT=$DATA/cohort-v0.4.19-changed-$ARM
[ -d "$COHORT" ] || COHORT=$DATA/cohort-v0.4.19-changed   # the first (fix2) run's cohort dir
OUT=${AB_OUT:-$LANE/$SUB}
export PYTHONPATH=$HERE/src

n=$(find "$COHORT" -name '*.xyz' | wc -l)
[ "$n" -gt 0 ] || { echo "ABORT: $COHORT is empty"; exit 1; }
for arm in off on; do
  d=$OUT/ab_$arm
  [ -f "$d/DONE" ] || { echo "ABORT: $d/DONE missing -- arm $arm still running (systemctl --user status oin-v0419-serab-$arm)"; exit 1; }
  r=$(find "$d/individual_reports" -name '*.json' | wc -l)
  [ "$r" = "$n" ] || { echo "ABORT: arm $arm has $r reports for a cohort of $n -- a shard died; re-run it with --continue UNDER THE SAME LEVERS"; exit 1; }
done

# OIN_H_FAITHFUL is a WRITER + READER lever (generation/metallogen_adapter.py preserves a bracketed
# [C] only when it is on), so the parse-back ruler -- which IS the reader -- must run under each
# arm's own configuration: the ON arm's strings read by the ON arm's reader, which is what a
# promotion ships. The same strings read by the SHIPPED reader are kept beside it
# (parseback_shippedreader.jsonl): that is what happens to a lever-written string in a build
# without the lever, and the report prints both.
cd "$HERE"
ALL_OFF="OIN_H_FAITHFUL=0 OIN_RC1_PROPAGATE=0 OIN_CAP_IGNORES_METAL=0 OIN_N_VALENCE_2=0 OIN_CANONICAL_RESONANCE=0 OIN_CANONICAL_CHARGES=0"
# The E2 arms stand on v0.4.19, where the three serializer levers ARE the shipped defaults: "off"
# there means only the E2 levers at 0 -- writing the promoted three to 0 would read the OFF arm's
# strings with the PRE-v0.4.19 reader and call the baseline something it is not.
case "$ARM" in e2*|res) ALL_OFF="OIN_N_VALENCE_2=0 OIN_CANONICAL_RESONANCE=0 OIN_CANONICAL_CHARGES=0";; esac
lev() { case "$1" in on) echo "$ON_LEVERS";; off) echo "$ALL_OFF";; esac; }
for arm in off on; do
  d=$OUT/ab_$arm
  [ -f "$d/bucket_report_honest.json" ] || $PY tools/roundtrip_bucket_report.py --results-dir "$d" --score honest > "$d/bucket_report_honest.log" 2>&1
  [ -f "$d/g_verdict.jsonl" ] || $PY tools/census/g_vs_input.py --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 8 > "$d/g_verdict.log" 2>&1
  if [ ! -f "$d/parseback.jsonl" ]; then
    env $(lev "$arm") $PY tools/census/string_sufficiency.py parseback --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 8 > "$d/parseback.log" 2>&1
  fi
  if [ ! -f "$d/parseback_gen.jsonl" ]; then
    env $(lev "$arm") $PY tools/census/string_sufficiency.py parseback --side gen --cohort "$COHORT" --sweep "$d" --out "$d" --cpu 8 > "$d/parseback_gen.log" 2>&1
  fi
  if [ "$arm" = on ] && [ ! -f "$d/parseback_shippedreader.jsonl" ]; then
    mkdir -p "$d/.shippedreader"
    env $(lev off) $PY tools/census/string_sufficiency.py parseback --cohort "$COHORT" --sweep "$d" --out "$d/.shippedreader" --cpu 8 > "$d/parseback_shippedreader.log" 2>&1
    mv "$d/.shippedreader/parseback.jsonl" "$d/parseback_shippedreader.jsonl"
  fi
done
sfx=""; [ "$ARM" = fix2 ] || sfx="_$ARM"
case "$ARM" in e2*|res) export OIN_AB_BASELINE=v0.4.19;; esac
$PY tools/v0419/serializer_ab_report.py all --ab "$OUT" --esc "$LANE/$ESC" --rescore "$LANE/$RESCORE" --single "$LANE/changed_single$sfx.jsonl" --out "$LANE/serializer_ab_report$sfx.json" | tee "$LANE/serializer_ab_report$sfx.txt"
