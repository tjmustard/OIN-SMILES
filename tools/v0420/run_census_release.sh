#!/usr/bin/env bash
# v0.4.20 -- the census on the RELEASE sweep (= the E2 candidate sweep,
# results-v0.4.19-e2-candidate-sweep: OIN_N_VALENCE_2 + OIN_CANONICAL_RESONANCE ON via the unit's
# environment at 86d9b0a4, promoted at 60288c14). Nothing else CPU-heavy may run beside it (load
# biases the probes' budgets).
#
#   tools/v0420/run_census_release.sh            # after post_candidate_sweep.sh; ~40 min on 10 CPUs
#
# tools/v0419/run_census_release.sh with the paths moved one release on. The encoder changed again,
# so C2 (e_selfconsistency) is the E2 lane's whole-cohort audit under the same two levers --
# results-v0.4.19-e2/e_selfconsistency_e2f.jsonl (5,000 x 11, 7dfaba98; src identical to the sweep's
# and to the promotion's). The CONTROL (the tool reproducing the frozen v0.4.17 census table) runs
# first as before. Movers are against the v0.4.19 release sweep.
set -euo pipefail

HERE=$(cd "$(dirname "$0")/../.." && pwd)
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
DATA=$MAIN/tmCAT-tmPHOTO_xyz_dataset
PY=$MAIN/.venv/bin/python
COHORT=$DATA/cohort-v0.4.5-5k
SWEEP=$DATA/results-v0.4.19-e2-candidate-sweep
OUT=$DATA/results-v0.4.20-release-census
OLD=$DATA/results-v0.4.19-release-census
CENSUS=$DATA/results-census
XF=$DATA/results-v0.4.19-e2
export PYTHONPATH=$HERE/src

[ -f "$SWEEP/POST_DONE" ] || { echo "ABORT: $SWEEP has not been post-processed (post_candidate_sweep.sh)"; exit 1; }
for f in g_verdict.jsonl parseback.jsonl bucket_report_honest.json; do
  [ -f "$SWEEP/$f" ] || { echo "ABORT: $SWEEP/$f missing"; exit 1; }
done
[ -z "$(git -C "$HERE" status --porcelain)" ] || { echo "REFUSING: $HERE is dirty"; exit 1; }
mkdir -p "$OUT/control"
git -C "$HERE" log -1 --format='%H %s' > "$OUT/COMMIT"
cd "$HERE"

step() { echo "== $1" | tee -a "$OUT/instruments.status"; }
step "control: the parameterised tool on the record sweep must reproduce the frozen census table"
$PY tools/census/attribution_table.py --write-to "$OUT/control" > "$OUT/control/attribution_control.txt" 2>&1
cmp <(zcat "$MAIN/measurements/v0.4.17-census/attribution_table.tsv.gz") "$OUT/control/attribution_table.tsv" \
  && echo "control: BYTE-IDENTICAL to the frozen census table" | tee -a "$OUT/instruments.status" \
  || { echo "ABORT: the control does not reproduce the frozen census table"; exit 1; }

step "attach"
$PY tools/attach_class_audit.py --results-dir "$SWEEP" --out "$OUT/attach_class_audit.json" > "$OUT/attach_class_audit.log" 2>&1; echo "rc=$?" >> "$OUT/instruments.status"
step "parseback gen"
$PY tools/census/string_sufficiency.py parseback --side gen --cohort "$COHORT" --sweep "$SWEEP" --out "$OUT" --cpu 10 > "$OUT/parseback_gen.log" 2>&1; echo "rc=$?" >> "$OUT/instruments.status"
step "pflags"
$PY tools/census/string_sufficiency.py pflags --charge-probe --cohort "$COHORT" --sweep "$SWEEP" --out "$OUT" --cpu 10 > "$OUT/pflags.log" 2>&1; echo "rc=$?" >> "$OUT/instruments.status"

step "attribution table"
$PY tools/census/attribution_table.py --sweep "$SWEEP" --out "$OUT" \
  --g-verdict "$SWEEP/g_verdict.jsonl" --parseback "$SWEEP/parseback.jsonl" \
  --parseback-gen "$OUT/parseback_gen.jsonl" --pflags "$OUT/pflags.jsonl" --attach "$OUT/attach_class_audit.json" \
  --e2 "$XF/e_selfconsistency_e2f.jsonl" \
  --mirror "$CENSUS/g_verdict_control-mirror.jsonl" --collisions "$CENSUS/collisions.json" \
  > "$OUT/attribution.txt" 2>&1; echo "rc=$?" >> "$OUT/instruments.status"

step "movers: structure differs from the v0.4.19 release sweep"
$PY - "$DATA" "$OUT" <<'PYEOF'
import hashlib, sys
from pathlib import Path
data, out = Path(sys.argv[1]), Path(sys.argv[2])
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
names = sorted(p.stem for p in (data / "cohort-v0.4.5-5k").glob("*.xyz"))
movers = [m for m in names if sha(data / "results-v0.4.19-candidate-sweep/structures" / f"{m}_generated.xyz") != sha(data / "results-v0.4.19-e2-candidate-sweep/structures" / f"{m}_generated.xyz")]
(out / "movers_vs_v0419.txt").write_text("\n".join(movers) + "\n")
print(f"movers {len(movers)} of {len(names)}")
PYEOF
step "reattribution diff vs the v0.4.19 release census"
$PY tools/census/reattribution_diff.py --old "$OLD/attribution_table.tsv" --new "$OUT/attribution_table.tsv" \
  --movers "$OUT/movers_vs_v0419.txt" --out "$OUT/reattribution_diff.json" > "$OUT/reattribution_diff.txt" 2>&1; echo "rc=$?" >> "$OUT/instruments.status"
echo ALLDONE >> "$OUT/instruments.status"
tail -60 "$OUT/attribution.txt"
