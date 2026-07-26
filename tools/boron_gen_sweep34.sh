#!/usr/bin/env bash
# Sweep GENERATION outcomes for the full 34-molecule boron-cage class (L4-boronfast, v0.4.7).
#
# tools/boron_gen_time.py already times ONE molecule and prints a JSON line; this wraps it in
# a HARD external `timeout` per molecule, because `embed_time_budget=self.timeout` is NOT a
# hard bound (docs/BORON_CAGE_v0.4.5.md SS10) -- it is checked only BETWEEN attempts in the
# embed loop, and a single attempt (an alt-cache MISS priming a fresh `option` via a full-complex
# PuLP/CBC solve) can itself run well past the requested budget. Only an external SIGKILL
# actually enforces a cap, so this script is that watchdog for the class-outcome sweep.
#
# Usage:
#   GEN_CAP=30 HARD_CAP=150 tools/boron_gen_sweep34.sh > tools/boron_gen_sweep34.jsonl
#
# Runs SERIALLY (1 worker), as required on a saturated box -- this is a class-outcome
# measurement (produced a structure: yes/no, exception type), not a timing benchmark, so
# wall-clock here is ADVISORY only.
set -uo pipefail

GEN_CAP="${GEN_CAP:-30}"
HARD_CAP="${HARD_CAP:-150}"
V="${OIN_PYTHON:-/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MOLS=(
  CAKBEW CAKBOG HAXJAS ICEZIC JABGAX PAQBOZ PAQCAM RANMUR RAWJEG RIWKAK RIWKEO ULODUU
  XUKRIF YIVLAQ GANYEZ HAXJOG JAFMIP JAFTAO JAFTES MAFSIY RONQET RONQOD YIBZIV COZCEZ
  GOHWOQ PAYTUH AVOFIB BEKLUA BEKMIP OZAREO MODZUA RANCIU RONPES RULBUV
)

export PYTHONPATH="$HERE/src"
export OIN_DATASET_DIR="${OIN_DATASET_DIR:-/home/tjmustard/Documents/GitHub/tmCat-tmPhoto/tmCAT-tmPHOTO_xyz_dataset}"

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
    # splice in the hard-timeout/wall fields without a JSON library dependency in the loop
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
