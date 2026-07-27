#!/usr/bin/env bash
# Control-group sweep: the 14 "silently wrong before, correctly cage-encoding now" molecules
# from docs/BORON_CAGE_v0.4.5.md SS5a (NOT the same population as the 34 encode_fail set --
# these 14 previously encoded a WRONG graph and passed the pipeline; OIN_BORON_CAGE=1 fixes
# their NOTATION but SS10's pipeline-arm gap was only sampled on the 34, not these 14).
#
# Purpose: the L4-boronfast fail-fast predicate flags on the B-B-B cage MOTIF alone, so it
# would also fire on any of these 14 if/when they hit generation. Before shipping it, check
# whether ANY of these 14 currently PRODUCE a 3D structure -- if even one does, the predicate
# must not be allowed to block it (task instruction: test the predicate against molecules that
# DO assemble, not only ones that don't).
set -uo pipefail

GEN_CAP="${GEN_CAP:-30}"
HARD_CAP="${HARD_CAP:-120}"
V="${OIN_PYTHON:-/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MOLS=(PEKQUU RAJNEY ULOFIK DUDTIG KIXXOF RAJNOI XIQKOY UYEJAK XIQLAL PEKQII VOFHUW CIDHAY SEMTOV VEJXOZ)

export PYTHONPATH="$HERE/src"
export OIN_DATASET_DIR="${OIN_DATASET_DIR:-/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset}"
export OIN_BORON_CAGE=1

for mol in "${MOLS[@]}"; do
  t0=$(date +%s.%N)
  out=$(GEN_CAP="$GEN_CAP" timeout -k 10 "$HARD_CAP" "$V" "$HERE/tools/boron_gen_time.py" "${mol}_comp_0" 2>/dev/null)
  rc=$?
  t1=$(date +%s.%N)
  wall=$(echo "$t1 - $t0" | bc)
  if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then
    printf '{"mol": "%s", "cap_s": %s, "hard_cap_s": %s, "hard_timeout": true, "wall_s": %.1f, "got_mol": false}\n' \
      "$mol" "$GEN_CAP" "$HARD_CAP" "$wall"
  elif [ -n "$out" ]; then
    echo "$out" | "$V" -c "
import json,sys
d=json.loads(sys.stdin.read())
d['hard_cap_s']=$HARD_CAP
d['hard_timeout']=False
d['wall_s']=round($wall,1)
print(json.dumps(d))
"
  else
    printf '{"mol": "%s", "cap_s": %s, "hard_cap_s": %s, "hard_timeout": true, "wall_s": %.1f, "got_mol": false, "note": "empty stdout, rc=%d"}\n' \
      "$mol" "$GEN_CAP" "$HARD_CAP" "$wall" "$rc"
  fi
done
