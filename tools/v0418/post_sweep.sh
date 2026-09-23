#!/usr/bin/env bash
# v0.4.18 -- read the finished sweep: honest buckets, the ruler on every generated structure,
# parse-back on every string, the two numbers, then the sweep against the A/B that predicted it.
#
#   tools/v0418/post_sweep.sh                                  # results-v0.4.18-sweep
#   SWEEP_TAG=v0.4.18-release SWEEP_RECORD=results-v0.4.18-sweep \
#     SWEEP_AB_ON=results-v0.4.18-eta-selection/ab_retarget SWEEP_PREDICTED=4305+45,3883+38 \
#     tools/v0418/post_sweep.sh                                # the release sweep, against the
#                                                              # sweep it replaces and the A/B
#
# Refuses an incomplete sweep. Completeness is REPORT COUNT == 5000: a shard that died early
# exits 0 like one that finished, and run_sweep.sh only warns.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
COHORT=$DATA/cohort-v0.4.5-5k
TAG=${SWEEP_TAG:-v0.4.18}
OUT=$DATA/results-$TAG-sweep
UNIT=oin-${TAG//./}-sweep
RECORD=$DATA/${SWEEP_RECORD:-results-v0.4.17-sweep}
AB_ON=$DATA/${SWEEP_AB_ON:-results-v0.4.18-eta-detached/ab_on}
PREDICTED=${SWEEP_PREDICTED:-4136+166,3730+150}
export PYTHONPATH=$HERE/src

[ -f "$OUT/DONE" ] || { echo "ABORT: $OUT/DONE missing -- still running? (systemctl --user status $UNIT)"; exit 1; }
r=$(find "$OUT/individual_reports" -name '*.json' | wc -l)
[ "$r" = 5000 ] || { echo "ABORT: $r reports, expected 5000 -- a shard died. Re-run it: tools/test_dataset_roundtrip.py --dataset-dir $COHORT --output-dir $OUT --shard I:6 --mol-timeout 300 --no-summary --continue (BLAS=1, NO OIN_* set), then re-write DONE"; exit 1; }
grep -q '"levers": {' "$OUT/run_config.json" && ! grep -q '"OIN_' "$OUT/run_config.json" || { echo "ABORT: run_config.json records an OIN_* lever -- this was not a shipped-defaults run"; exit 1; }

cd "$HERE"
$PY tools/v0417/sweep_two_numbers.py --control
$PY tools/roundtrip_bucket_report.py --results-dir "$OUT" --score honest > "$OUT/bucket_report_honest.log" 2>&1
$PY tools/census/g_vs_input.py --cohort "$COHORT" --sweep "$OUT" --out "$OUT" --cpu 10 > "$OUT/g_verdict.log" 2>&1
$PY tools/census/string_sufficiency.py parseback --cohort "$COHORT" --sweep "$OUT" --out "$OUT" --cpu 10 > "$OUT/parseback.log" 2>&1
# NOTE the "(record ...)" labels this prints are the v0.4.14 record's; the comparison that matters
# here -- against results-v0.4.17-sweep and against the A/B -- is the next tool's.
$PY tools/v0417/sweep_two_numbers.py --sweep "$OUT" | tee "$OUT/two_numbers.txt"
$PY tools/v0418/sweep_vs_ab.py --sweep "$OUT" --record "$RECORD" --ab-on "$AB_ON" --predicted "$PREDICTED" --out "$OUT/sweep_vs_ab.json" | tee "$OUT/sweep_vs_ab.txt"
