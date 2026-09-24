#!/usr/bin/env bash
# v0.4.18 -- the full sweep the eta construction promotion owes, launched detached.
#
# A thin wrapper, on purpose: the sweep itself is tools/run_sweep.sh, the project's own
# version-controlled launcher (interpreter resolved and rdkit pin REQUIRED, PYTHONPATH pinned,
# run_config.json written). This only supplies what that script does not:
#
#   * a systemd --user SERVICE unit, so the run survives the session (`--scope` rejects
#     OOMPolicy; a service accepts it) and an OOM kill takes one child, not the sweep;
#   * the sweep of record's thread caps (results-v0.4.14-sweep/RUN.md): BLAS=1, six shards on
#     twelve cores. OIN3DGenerator(timeout=) is ADVISORY, so CPU starvation shrinks the pool and
#     UNDERSTATES byte_exact -- a contended sweep is a biased sweep;
#   * NO OIN_* VARIABLE AT ALL. The point of this run is the SHIPPED DEFAULT. A --user unit does
#     not inherit the login shell's environment, and nothing below passes a lever, so the only
#     way the eta levers are on is that levers._DEFAULT_ON says so. run_config.json's "levers"
#     block must come back EMPTY -- check it.
#
# Like-for-like with the sweep of record: same cohort, 6 shards (1-based), --mol-timeout 300.
#
#   tools/v0418/launch_sweep.sh           # refuses a dirty tree: every report stamps HEAD
#   SWEEP_TAG=v0.4.18-release tools/v0418/launch_sweep.sh
#                                         # -> results-v0.4.18-release-sweep, unit oin-v0418-release-sweep
#                                         #    (the sweep OIN_ETA_RETARGET's promotion owes)
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
DATA=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset
COHORT=$DATA/cohort-v0.4.5-5k
TAG=${SWEEP_TAG:-v0.4.18}
OUT=$DATA/results-$TAG-sweep
UNIT=oin-${TAG//./}-sweep

[ -z "$(git -C "$HERE" status --porcelain)" ] || { echo "ABORT: $HERE is dirty -- the reports would stamp a commit that is not the code that ran"; exit 1; }
[ "$(find "$COHORT" -xtype l | wc -l)" = 0 ] || { echo "ABORT: dangling cohort symlinks -- restore the dataset first"; exit 1; }
[ "$(find "$COHORT" -name '*.xyz' | wc -l)" = 5000 ] || { echo "ABORT: cohort is not 5000 molecules"; exit 1; }
[ ! -e "$OUT/individual_reports" ] || { echo "ABORT: $OUT already has reports -- a second launch would mix two runs"; exit 1; }
if env | grep -q '^OIN_'; then echo "note: this shell has OIN_* set; the unit will NOT inherit them:"; env | grep '^OIN_'; fi

mkdir -p "$OUT"
systemd-run --user --unit="$UNIT" \
  --description="$TAG full sweep, shipped defaults (levers._DEFAULT_ON as committed)" \
  -p OOMPolicy=continue -p MemoryMax=14G --working-directory="$HERE" \
  -E PATH="$PATH" -E HOME="$HOME" \
  -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E NUMEXPR_NUM_THREADS=1 \
  /bin/bash -c "tools/run_sweep.sh $COHORT $OUT 6 300 > $OUT/launch.log 2>&1; echo \"#DONE \$(ls $OUT/individual_reports | wc -l)\" > $OUT/DONE"
echo "launched $UNIT from $(git -C "$HERE" rev-parse --short HEAD). Finished when $OUT/DONE exists; then SWEEP_TAG=$TAG tools/v0418/post_sweep.sh"
