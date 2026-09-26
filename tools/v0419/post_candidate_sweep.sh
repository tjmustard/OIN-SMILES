#!/usr/bin/env bash
# v0.4.19 -- read the finished CANDIDATE sweep (the three serializer levers ON via the unit's
# environment): honest buckets, the neutral ruler on every generated structure, parse-back on
# every string UNDER THE LEVERS, the two numbers, then the sweep against the projection that
# predicted it, row by row.
#
#   SWEEP_TAG=v0.4.19-candidate tools/v0419/post_candidate_sweep.sh
#   SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS="OIN_N_VALENCE_2=1 OIN_CANONICAL_RESONANCE=1" \
#     SWEEP_SHIPPED="OIN_N_VALENCE_2=0 OIN_CANONICAL_RESONANCE=0" SWEEP_PROJECTION=e2f \
#     tools/v0419/post_candidate_sweep.sh
#   (the E2 lane stands ON v0.4.19: its "shipped reader" leaves the promoted three at their default)
#
# Differences from tools/v0418/post_sweep.sh, each on purpose:
#   * run_config.json's "levers" block must be EXACTLY the three (post_sweep.sh requires it empty);
#   * OIN_H_FAITHFUL is a WRITER + READER lever (generation/metallogen_adapter.py keeps a bracketed
#     [C] only while it is on), so the parse-back ruler -- which IS the reader -- runs under the
#     three levers: that is what a promotion ships. The same strings read by the SHIPPED reader are
#     kept beside it (parseback_shippedreader.jsonl) so the coupling is visible, never hidden;
#   * the comparison is tools/v0419/candidate_vs_projection.py against the fix3 projection
#     (release sweep + offline re-score + the live changed-set arm), not sweep_vs_ab.py's eta logic.
# Refuses an incomplete sweep: completeness is REPORT COUNT == 5000.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
COHORT=$DATA/cohort-v0.4.5-5k
TAG=${SWEEP_TAG:-v0.4.19-candidate}
OUT=$DATA/results-$TAG-sweep
LEVERS=${SWEEP_LEVERS:-"OIN_H_FAITHFUL=1 OIN_RC1_PROPAGATE=1 OIN_CAP_IGNORES_METAL=1"}
SHIPPED=${SWEEP_SHIPPED:-"OIN_H_FAITHFUL=0 OIN_RC1_PROPAGATE=0 OIN_CAP_IGNORES_METAL=0"}
PROJECTION=${SWEEP_PROJECTION:-fix3}
export PYTHONPATH=$HERE/src

[ -f "$OUT/DONE" ] || { echo "ABORT: $OUT/DONE missing -- still running?"; exit 1; }
r=$(find "$OUT/individual_reports" -name '*.json' | wc -l)
[ "$r" = 5000 ] || { echo "ABORT: $r reports, expected 5000 -- a shard died"; exit 1; }
nlev=0
for kv in $LEVERS; do
  k=${kv%%=*}; nlev=$((nlev + 1))
  grep -q "\"$k\": \"1\"" "$OUT/run_config.json" || { echo "ABORT: run_config.json does not record $k=1"; exit 1; }
done
[ "$(grep -c '"OIN_' "$OUT/run_config.json")" = "$nlev" ] || { echo "ABORT: run_config.json records a lever beyond $LEVERS"; exit 1; }

cd "$HERE"
$PY tools/v0417/sweep_two_numbers.py --control
[ -f "$OUT/bucket_report_honest.json" ] || $PY tools/roundtrip_bucket_report.py --results-dir "$OUT" --score honest > "$OUT/bucket_report_honest.log" 2>&1
[ -f "$OUT/g_verdict.jsonl" ] || $PY tools/census/g_vs_input.py --cohort "$COHORT" --sweep "$OUT" --out "$OUT" --cpu 10 > "$OUT/g_verdict.log" 2>&1
[ -f "$OUT/parseback.jsonl" ] || env $LEVERS $PY tools/census/string_sufficiency.py parseback --cohort "$COHORT" --sweep "$OUT" --out "$OUT" --cpu 10 > "$OUT/parseback.log" 2>&1
[ -f "$OUT/parseback_gen.jsonl" ] || env $LEVERS $PY tools/census/string_sufficiency.py parseback --side gen --cohort "$COHORT" --sweep "$OUT" --out "$OUT" --cpu 10 > "$OUT/parseback_gen.log" 2>&1
if [ ! -f "$OUT/parseback_shippedreader.jsonl" ]; then
  mkdir -p "$OUT/.shippedreader"
  env $SHIPPED $PY tools/census/string_sufficiency.py parseback --cohort "$COHORT" --sweep "$OUT" --out "$OUT/.shippedreader" --cpu 10 > "$OUT/parseback_shippedreader.log" 2>&1
  mv "$OUT/.shippedreader/parseback.jsonl" "$OUT/parseback_shippedreader.jsonl"
fi
$PY tools/v0417/sweep_two_numbers.py --sweep "$OUT" | tee "$OUT/two_numbers.txt"
$PY tools/v0419/candidate_vs_projection.py --sweep "$OUT" --projection "$PROJECTION" --out "$OUT/candidate_vs_projection.json" | tee "$OUT/candidate_vs_projection.txt"
echo "#DONE post" > "$OUT/POST_DONE"
