#!/bin/bash
# v0.4.7 byte-identity gate. Two independent arms, both diffing SHA256 of the OIN
# *string* against a frozen golden manifest -- never the generated XYZ (see
# tools/gate_arm2_roundtrip_one.py's docstring for why that would be the wrong,
# too-strong gate object).
#
#   tools/gate_v047.sh arm1 [--fixtures-dir DIR] [--golden PATH] [--out PATH]
#   tools/gate_v047.sh arm2 [--cohort-dir DIR] [--golden PATH] [--out PATH] [--timeout S]
#   tools/gate_v047.sh both [same options as above]
#
# ARM 1 (encoder, 61 fixtures): tools/gate_arm1_encode.py, one process, memo-cleared
#   between molecules (see that script's docstring). Fast -- runs to completion here.
#
# ARM 2 (round trip, the frozen slow-cohort): tools/gate_arm2_roundtrip_one.py, ONE
#   SUBPROCESS PER MOLECULE (see that script's docstring for why: several OIN_* levers
#   and module-level caches are frozen at import time, and a fresh interpreter per
#   molecule is the only isolation guarantee that does not depend on enumerating every
#   such cache correctly). This arm is EXPENSIVE by construction -- the cohort is
#   exactly the molecules known to be slow -- so do not run it against the full
#   cohort on a loaded box without knowing what you are signing up for.
#
# MANDATORY DISCIPLINE (all learned the hard way on this project -- see MEMORY.md):
#   * Every result line is written with a real `>>` file append (not buffered inside
#     a long-lived process) and the per-molecule python scripts flush(=True) every
#     print, so a `timeout`/kill mid-run cannot discard output that looks like
#     agreement once `sort`ed.
#   * The `#DONE <n>` sentinel is REQUIRED and its denominator is checked BEFORE any
#     comparison is trusted -- an empty or truncated results file must never look
#     like consensus.
#   * Never `2>/dev/null` -- stderr is always let through to the log file.
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Resolution order matches tools/run_sweep.sh: explicit override, this repo's own
# venv, then a sibling checkout's. Worktrees have no .venv of their own -- project
# rule is to use the MAIN checkout's pinned venv (rdkit 2025.9.3) and never `uv sync`.
if [ -n "${OIN_SWEEP_PYTHON:-}" ]; then
    PY="$OIN_SWEEP_PYTHON"
elif [ -x "$REPO/.venv/bin/python" ]; then
    PY="$REPO/.venv/bin/python"
else
    for cand in "$(dirname "$REPO")"/*/.venv/bin/python; do
        [ -x "$cand" ] && PY="$cand" && break
    done
fi
if [ -z "${PY:-}" ] || [ ! -x "$PY" ]; then
    echo "error: no python found. Set OIN_SWEEP_PYTHON to the MAIN checkout's" >&2
    echo "       .venv/bin/python (rdkit is pinned there; never uv sync in a worktree)." >&2
    exit 1
fi
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

ARM="${1:-}"
shift || true

FIXTURES_DIR="$REPO/tests/fixtures"
ARM1_GOLDEN="$REPO/tools/gate_v047_arm1_golden.tsv"
COHORT_DIR=""
ARM2_GOLDEN="$REPO/tools/gate_v047_arm2_golden.tsv"
GEN_TIMEOUT=300
OUT=""

while [ $# -gt 0 ]; do
    case "$1" in
        --fixtures-dir) FIXTURES_DIR="$2"; shift 2 ;;
        --cohort-dir) COHORT_DIR="$2"; shift 2 ;;
        --golden)
            if [ "$ARM" = "arm2" ]; then ARM2_GOLDEN="$2"; else ARM1_GOLDEN="$2"; fi
            shift 2 ;;
        --timeout) GEN_TIMEOUT="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        *) echo "error: unrecognized option $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
run_arm1() {
    local out="${OUT:-$REPO/spec/handoffs/v0.4.7/gate_arm1_result.tsv}"
    mkdir -p "$(dirname "$out")"
    echo "[gate/arm1] fixtures_dir=$FIXTURES_DIR out=$out" >&2

    "$PY" "$REPO/tools/gate_arm1_encode.py" --fixtures-dir "$FIXTURES_DIR" > "$out"
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "[gate/arm1] FAIL: gate_arm1_encode.py exited $rc" >&2
        return 1
    fi

    local done_line n_done
    done_line="$(grep -E '^#DONE ' "$out" | tail -1)"
    if [ -z "$done_line" ]; then
        echo "[gate/arm1] FAIL: no #DONE sentinel in $out -- run was truncated/killed" >&2
        return 1
    fi
    n_done="$(echo "$done_line" | awk '{print $2}')"
    if [ "$n_done" -ne 61 ]; then
        echo "[gate/arm1] FAIL: #DONE $n_done, expected 61 -- short run, refusing to trust it" >&2
        return 1
    fi
    echo "[gate/arm1] sentinel OK: #DONE $n_done" >&2

    local manifest_line
    manifest_line="$(grep -E '^# MANIFEST_SHA256=' "$out" | tail -1)"
    echo "[gate/arm1] $manifest_line" >&2

    if [ ! -f "$ARM1_GOLDEN" ]; then
        echo "[gate/arm1] no golden at $ARM1_GOLDEN -- nothing to compare against (first run?)" >&2
        return 0
    fi

    local golden_manifest fresh_manifest
    golden_manifest="$(grep -E '^# MANIFEST_SHA256=' "$ARM1_GOLDEN" | tail -1)"
    fresh_manifest="$manifest_line"
    if [ "$golden_manifest" != "$fresh_manifest" ]; then
        echo "[gate/arm1] MISMATCH" >&2
        echo "[gate/arm1]   golden: $golden_manifest" >&2
        echo "[gate/arm1]   fresh:  $fresh_manifest" >&2
        echo "[gate/arm1] per-line diff (data rows only, sorted):" >&2
        diff <(grep -vE '^#' "$ARM1_GOLDEN" | sort) <(grep -vE '^#' "$out" | sort) >&2
        return 1
    fi
    echo "[gate/arm1] PASS -- byte-identical to golden" >&2
    return 0
}

# ---------------------------------------------------------------------------
run_arm2() {
    if [ -z "$COHORT_DIR" ]; then
        echo "error: arm2 requires --cohort-dir" >&2
        return 1
    fi
    if [ ! -d "$COHORT_DIR" ]; then
        echo "[gate/arm2] FAIL: cohort dir $COHORT_DIR not found" >&2
        return 1
    fi

    local out="${OUT:-$REPO/spec/handoffs/v0.4.7/gate_arm2_result.tsv}"
    mkdir -p "$(dirname "$out")"
    : > "$out"  # truncate/create -- real file, one `>>` append per molecule below

    local names
    names="$(find "$COHORT_DIR" -maxdepth 1 -name '*.xyz' -printf '%f\n' | sed 's/\.xyz$//' | sort)"
    local expected
    expected="$(echo "$names" | grep -c .)"
    if [ "$expected" -eq 0 ]; then
        echo "[gate/arm2] FAIL: 0 molecules found under $COHORT_DIR -- refusing an empty corpus" >&2
        return 1
    fi
    echo "[gate/arm2] cohort_dir=$COHORT_DIR molecules=$expected timeout=${GEN_TIMEOUT}s out=$out" >&2
    echo "[gate/arm2] one subprocess per molecule -- this arm is EXPENSIVE by construction" >&2

    local i=0
    while IFS= read -r name; do
        [ -z "$name" ] && continue
        i=$((i + 1))
        echo "[gate/arm2] ($i/$expected) $name ..." >&2
        "$PY" "$REPO/tools/gate_arm2_roundtrip_one.py" \
            --cohort-dir "$COHORT_DIR" --molecule "$name" --timeout "$GEN_TIMEOUT" \
            >> "$out"
    done <<< "$names"

    sort -o "$out" "$out"
    echo "#DONE $i" >> "$out"

    local n_done
    n_done="$(grep -E '^#DONE ' "$out" | tail -1 | awk '{print $2}')"
    if [ "$n_done" -ne "$expected" ]; then
        echo "[gate/arm2] FAIL: #DONE $n_done, expected $expected -- short run, refusing to trust it" >&2
        return 1
    fi
    echo "[gate/arm2] sentinel OK: #DONE $n_done" >&2

    if [ ! -f "$ARM2_GOLDEN" ]; then
        echo "[gate/arm2] no golden at $ARM2_GOLDEN -- nothing to compare against (first run?)" >&2
        return 0
    fi

    # Compare per-molecule sha256(smiles_1)/sha256(smiles_2) columns (2 and 3) keyed
    # by name (column 1). The gate object is the STRING hash, never xyz_sha256
    # (column 7, an observation-only column -- see gate_arm2_roundtrip_one.py).
    local mismatches=0 missing=0
    while IFS=$'\t' read -r g_name g_sha1 g_sha2 _rest; do
        [ -z "$g_name" ] && continue
        local row
        row="$(grep -E "^${g_name}"$'\t' "$out" || true)"
        if [ -z "$row" ]; then
            echo "[gate/arm2] MISSING: $g_name has a golden row but no fresh result" >&2
            missing=$((missing + 1))
            continue
        fi
        local f_sha1 f_sha2
        f_sha1="$(echo "$row" | cut -f2)"
        f_sha2="$(echo "$row" | cut -f3)"
        if [ "$f_sha1" != "$g_sha1" ] || [ "$f_sha2" != "$g_sha2" ]; then
            echo "[gate/arm2] MISMATCH: $g_name" >&2
            echo "[gate/arm2]   golden sha1=$g_sha1 sha2=$g_sha2" >&2
            echo "[gate/arm2]   fresh  sha1=$f_sha1 sha2=$f_sha2" >&2
            echo "[gate/arm2]   fresh row: $row" >&2
            mismatches=$((mismatches + 1))
        fi
    done < <(grep -vE '^#' "$ARM2_GOLDEN")

    if [ $mismatches -gt 0 ] || [ $missing -gt 0 ]; then
        echo "[gate/arm2] FAIL: $mismatches mismatch(es), $missing missing" >&2
        return 1
    fi
    echo "[gate/arm2] PASS -- all $expected molecules byte-identical to golden" >&2
    return 0
}

# ---------------------------------------------------------------------------
case "$ARM" in
    arm1) run_arm1 ;;
    arm2) run_arm2 ;;
    both)
        run_arm1
        rc1=$?
        run_arm2
        rc2=$?
        [ $rc1 -eq 0 ] && [ $rc2 -eq 0 ]
        ;;
    *)
        echo "usage: $0 <arm1|arm2|both> [--fixtures-dir DIR] [--cohort-dir DIR] [--golden PATH] [--timeout S] [--out PATH]" >&2
        exit 1
        ;;
esac
