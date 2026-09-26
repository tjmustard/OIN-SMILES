#!/usr/bin/env bash
# v0.4.19 -- the CANDIDATE sweep: the full 5,000 through the real harness with the three measured
# levers ON, so the owner's promotion decision rests on a sweep and not on a projection.
#
#   OIN_H_FAITHFUL=1  OIN_RC1_PROPAGATE=1  OIN_CAP_IGNORES_METAL=1
#
# Everything else is tools/v0418/launch_sweep.sh verbatim (a --user SERVICE unit, BLAS=1, six
# 1-based shards, --mol-timeout 300, same cohort) -- that launcher passes NO OIN_* variable on
# purpose (a --user unit does not inherit the login shell), so this one exists only to add the
# three -E flags. tools/run_sweep.sh copies every OIN_* it sees into run_config.json's "levers"
# block: that block must come back with exactly these three -- check it before trusting a number.
#
# The prediction this sweep tests (docs/agentic-notes/v0.4.19/SERIALIZER_LANE.md section 6b, from
# the release sweep + offline re-score + the three live changed-set arms): 4,378 / 3,951 of 5,000 =
# 87.56% self-consistent / 79.02% VERIFIED. Budget-boundary rows will differ; the levers are the
# only intended difference.
#
#   tools/v0419/launch_candidate_sweep.sh     # -> results-v0.4.19-candidate-sweep, unit oin-v0419-candidate-sweep
#   then: SWEEP_TAG=v0.4.19-candidate tools/v0419/post_candidate_sweep.sh
#
# The E2 lane's candidate (v0.4.19 defaults + the two perception levers; projection =
# results-v0.4.19-e2/serializer_ab_report_e2f.json, 4,466 / 4,025 = 89.32% / 80.50%):
#   SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS="-E OIN_N_VALENCE_2=1 -E OIN_CANONICAL_RESONANCE=1" \
#       tools/v0419/launch_candidate_sweep.sh
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
DATA=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset
COHORT=$DATA/cohort-v0.4.5-5k
TAG=${SWEEP_TAG:-v0.4.19-candidate}
OUT=$DATA/results-$TAG-sweep
UNIT=oin-${TAG//./}-sweep
LEVERS=${SWEEP_LEVERS:-"-E OIN_H_FAITHFUL=1 -E OIN_RC1_PROPAGATE=1 -E OIN_CAP_IGNORES_METAL=1"}

[ -z "$(git -C "$HERE" status --porcelain)" ] || { echo "ABORT: $HERE is dirty -- the reports would stamp a commit that is not the code that ran"; exit 1; }
if systemctl --user show-environment | grep -q '^OIN_'; then echo "ABORT: the user manager exports an OIN_* variable; the unit would see more than $LEVERS"; exit 1; fi
[ "$(find "$COHORT" -xtype l | wc -l)" = 0 ] || { echo "ABORT: dangling cohort symlinks -- restore the dataset first"; exit 1; }
[ "$(find "$COHORT" -name '*.xyz' | wc -l)" = 5000 ] || { echo "ABORT: cohort is not 5000 molecules"; exit 1; }
[ ! -e "$OUT/individual_reports" ] || { echo "ABORT: $OUT already has reports -- a second launch would mix two runs"; exit 1; }
if systemctl --user list-units --no-legend 'oin-*' | grep -q .; then echo "ABORT: another oin-* unit is running; load biases accuracy"; systemctl --user list-units --no-legend 'oin-*'; exit 1; fi

mkdir -p "$OUT"
systemd-run --user --unit="$UNIT" \
  --description="$TAG full sweep, levers ON ($LEVERS)" \
  -p OOMPolicy=continue -p MemoryMax=14G --working-directory="$HERE" \
  -E PATH="$PATH" -E HOME="$HOME" \
  $LEVERS \
  -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E NUMEXPR_NUM_THREADS=1 \
  /bin/bash -c "tools/run_sweep.sh $COHORT $OUT 6 300 > $OUT/launch.log 2>&1; echo \"#DONE \$(ls $OUT/individual_reports | wc -l)\" > $OUT/DONE"
echo "launched $UNIT from $(git -C "$HERE" rev-parse --short HEAD). Finished when $OUT/DONE exists; then SWEEP_TAG=$TAG tools/v0419/post_candidate_sweep.sh (same SWEEP_* variables)"
