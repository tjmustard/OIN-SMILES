#!/usr/bin/env bash
# Staged runner for the CREST conformer-invariance sweep
# (tools/conformer_invariance_crest.py). Handles crest discovery, sane defaults, a
# persistent (non-/tmp) output dir, and an optional systemd-run wrapper with OOM guards
# for the long full sweep.
#
# Usage:
#   bash tools/run_conformer_crest_sweep.sh smoke            # 3 tiny structures, foreground
#   bash tools/run_conformer_crest_sweep.sh small            # 8 smallest structures, foreground
#   bash tools/run_conformer_crest_sweep.sh full             # all 30, foreground
#   bash tools/run_conformer_crest_sweep.sh full --systemd   # all 30 as a transient user
#                                                            # service (survives shell/harness
#                                                            # exit; OOMPolicy=continue)
#
# Tunables via env:
#   OIN_VENV_PY       python that imports oinsmiles (default: main-checkout .venv)
#   OIN_CREST_THREADS crest -T value                (default 4)
#   OIN_CREST_TIMEOUT per-structure seconds         (default 1800)
#   OIN_CREST_METHOD  crest method flag             (default gfnff)
#   OIN_CREST_MAXCONF conformers encoded/structure  (default 10)
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
VENV_PY="${OIN_VENV_PY:-/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python}"
THREADS="${OIN_CREST_THREADS:-4}"
TIMEOUT="${OIN_CREST_TIMEOUT:-1800}"
METHOD="${OIN_CREST_METHOD:-gfnff}"
MAXCONF="${OIN_CREST_MAXCONF:-10}"
PREOPT="${OIN_CREST_PREOPT:-none}"          # none|gxtb|gfn2|gfnff (xtb pre-opt before CREST)
REOPT="${OIN_CREST_REOPT:-none}"            # none|gxtb|gfn2|gfnff (xtb re-opt of EACH conformer)
CREST_ARGS="${OIN_CREST_ARGS:-}"            # e.g. "--noreftopo --mquick"
XTB_BIN="${OIN_XTB_BIN:-/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/xtb}"
KEEPTMP="${OIN_CREST_KEEPTMP:-0}"           # 1 = keep per-structure CREST/pre-opt scratch
OUT="${OIN_CREST_OUT:-$REPO/conformer_crest_sweep}"  # persistent + gitignored; never /tmp

STAGE="${1:-}"; shift || true
case "$STAGE" in
    smoke) SCOPE=(--only CisPlatin,TransPlatin,XIZFAQ_comp_0 --max-confs 5) ;;
    small) SCOPE=(--limit 8 --max-confs "$MAXCONF") ;;
    full)  SCOPE=(--max-confs "$MAXCONF") ;;
    *) echo "usage: $0 {smoke|small|full} [--systemd]" >&2; exit 2 ;;
esac

# --- Discover crest: PATH -> conda-forge env 'crest' -> repo-local .crest/bin --------
if ! command -v crest >/dev/null 2>&1; then
    for MGR in micromamba mamba conda; do
        command -v "$MGR" >/dev/null 2>&1 || continue
        CB="$("$MGR" run -n crest bash -c 'command -v crest' 2>/dev/null || true)"
        if [ -n "$CB" ]; then export PATH="$(dirname "$CB"):$PATH"; break; fi
    done
fi
if ! command -v crest >/dev/null 2>&1 && [ -x "$REPO/.crest/bin/crest" ]; then
    export PATH="$REPO/.crest/bin:$PATH"
fi
if ! command -v crest >/dev/null 2>&1; then
    # CREST is an OPTIONAL external tool, not a dependency. A missing binary is not an
    # error: skip cleanly (exit 0) so this runner is CI-safe and hermetic when CREST is
    # absent, mirroring the shutil.which gate inside conformer_invariance_crest.py.
    echo "CREST is not installed (no 'crest' on PATH, no conda-forge 'crest' env, no"
    echo "$REPO/.crest/bin/crest). This conformer cross-check is OPTIONAL and not a"
    echo "dependency of OIN-SMILES -- skipping cleanly. Install it with"
    echo "'bash tools/install_crest.sh' to enable it, then re-run."
    exit 0
fi
echo "crest: $(command -v crest)  ($(crest --version 2>/dev/null | grep -i version | tr -s ' '))"

mkdir -p "$OUT"
CMD=("$VENV_PY" "$REPO/tools/conformer_invariance_crest.py"
     "${SCOPE[@]}" --method "$METHOD" --threads "$THREADS"
     --crest-timeout "$TIMEOUT" --preopt "$PREOPT" --reopt "$REOPT"
     --xtb-bin "$XTB_BIN" --out-dir "$OUT")
[ -n "$CREST_ARGS" ] && CMD+=(--crest-args "$CREST_ARGS")
[ "$KEEPTMP" = "1" ] && CMD+=(--keep-tmp)
echo "config: preopt=$PREOPT reopt=$REOPT method=$METHOD crest-args='${CREST_ARGS:-none}' keep-tmp=$KEEPTMP timeout=${TIMEOUT}s out=$OUT"

if [ "${1:-}" = "--systemd" ]; then
    # Long full sweep: transient user service, survives shell/harness exit. OOMPolicy=continue
    # so one runaway forked child cannot take down the whole run (lesson from prior sweeps).
    UNIT="${OIN_CREST_UNIT:-crest-conformer-sweep}"
    systemctl --user reset-failed "${UNIT}.service" 2>/dev/null || true
    echo "Launching '$UNIT' via systemd-run --user (OOMPolicy=continue, MemoryMax=14G)."
    echo "Follow: journalctl --user -u ${UNIT}.service -f    Output: $OUT"
    exec systemd-run --user --unit="$UNIT" \
        -p OOMPolicy=continue -p MemoryMax=14G \
        --working-directory="$REPO" \
        env "PATH=$PATH" "${CMD[@]}"
else
    echo "Running foreground. Output: $OUT"
    exec "${CMD[@]}"
fi
