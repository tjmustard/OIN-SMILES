#!/bin/bash
# v0.4.7 byte-identity gate. Two independent arms, both diffing SHA256 of the OIN
# *string* against a frozen golden manifest -- never the generated XYZ (see
# tools/gate_arm2_roundtrip_one.py's docstring for why that would be the wrong,
# too-strong gate object).
#
#   tools/gate_v047.sh arm1 [--fixtures-dir DIR] [--golden PATH] [--out PATH]
#   tools/gate_v047.sh arm2 [--cohort-dir DIR] [--golden PATH] [--out PATH] [--timeout S]
#                           [--hard-timeout S] [--shard i:n] [--band NAME]
#   tools/gate_v047.sh both [same options as above]
#
# v0.4.9 additions, all on ARM 2 and all forced by the stratified runtime cohort
# (tools/gate_v049_arm2_golden.tsv, 325 rows, ~10 CPU-h serially):
#   --hard-timeout  a REAL SIGKILL per molecule. `--timeout` is only the generator's
#                   advisory budget, so before this the gate itself was unbounded.
#   --shard i:n     ONE-BASED, matching test_dataset_roundtrip.py. 10 CPU-h serial is not
#                   a thing anyone runs twice for a reproducibility check.
#   --band NAME     one runtime stratum, so the fast bands can be gated routinely.
# The golden may also carry NO_STRUCTURE@Ns / NO_ENCODE@Ns sentinels in place of a hash --
# see the SENTINELS comment in run_arm2 for why those must not be compared as hashes.
#
# ARM 1 (encoder, 62 fixtures): tools/gate_arm1_encode.py, one process, memo-cleared
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
# ⚠ THE INTERPRETER IS RESOLVED, NEVER GUESSED.
#
# The previous fallback globbed `$(dirname $REPO)/*/.venv/bin/python` because worktrees have
# no .venv of their own. That glob matches every project sharing the parent directory. Run
# from a v0.4.9 worktree under ~/Documents/GitHub it selected `EtaCatalysis/.venv` -- an
# unrelated project with no rdkit at all -- and then `EtaTMCSMILES/.venv`, which has rdkit
# **2025.09.2** against this project's pinned **2025.9.3**. The first failure is loud but
# misdiagnosed (a short `#DONE`, which reads as "the run was truncated"); the second is
# worse, because a gate run on a different rdkit produces MISMATCHes that look like code
# regressions. Cross-worktree rdkit drift has already cost this project an A/B.
#
# A worktree resolves to its own main checkout deterministically via `--git-common-dir`, so
# there is nothing to guess. If that does not yield a usable interpreter, this fails.
_rdkit_version() { "$1" -c 'import rdkit; print(rdkit.__version__)' 2>/dev/null; }

PIN="$(sed -n 's/.*"rdkit==\([^"]*\)".*/\1/p' "$REPO/pyproject.toml" | head -1)"

MAIN_CHECKOUT="$(dirname "$(git -C "$REPO" rev-parse --git-common-dir 2>/dev/null || echo /nonexistent/.git)")"
PY=""
_tried=""
for cand in "${OIN_SWEEP_PYTHON:-}" "$REPO/.venv/bin/python" "$MAIN_CHECKOUT/.venv/bin/python"; do
    [ -z "$cand" ] || [ ! -x "$cand" ] && { _tried="$_tried ${cand:-<unset>}"; continue; }
    _tried="$_tried $cand"
    [ -n "$(_rdkit_version "$cand")" ] && { PY="$cand"; break; }
done
if [ -z "$PY" ]; then
    echo "error: no python with an importable rdkit found. Set OIN_SWEEP_PYTHON to the MAIN" >&2
    echo "       checkout's .venv/bin/python (rdkit is pinned ==$PIN there; never 'uv sync'" >&2
    echo "       in a worktree). Tried:$_tried" >&2
    exit 1
fi
RDKIT_V="$(_rdkit_version "$PY")"
echo "[gate] python=$PY (rdkit $RDKIT_V)" >&2
# Normalize 2025.09.2 vs 2025.9.2 before comparing -- zero-padding differs between the
# pyproject pin and rdkit.__version__ and a spurious warning trains people to ignore it.
_norm() { echo "$1" | awk -F. '{for(i=1;i<=NF;i++){printf "%s%d", (i>1?".":""), $i}; print ""}'; }
if [ -n "$PIN" ] && [ "$(_norm "$RDKIT_V")" != "$(_norm "$PIN")" ]; then
    echo "error: rdkit $RDKIT_V != pinned $PIN. A byte-identity gate on a different rdkit" >&2
    echo "       reports MISMATCHes that are version drift, not code changes. Refusing." >&2
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
# Hard wall-clock kill per molecule. `--timeout` is only the generator's ADVISORY budget
# (OIN3DGenerator(timeout=) is checked between embed attempts, never inside one -- that is
# v0.4.9 Lane 1's whole subject), so without this the gate itself is unbounded. Default is
# 1.5x the generator budget: enough headroom that a molecule finishing legitimately late is
# not cut, small enough that the arm has a stated worst case. 0 disables.
HARD_TIMEOUT=""
SHARD=""
BAND=""
OUT=""

while [ $# -gt 0 ]; do
    case "$1" in
        --fixtures-dir) FIXTURES_DIR="$2"; shift 2 ;;
        --cohort-dir) COHORT_DIR="$2"; shift 2 ;;
        --golden)
            if [ "$ARM" = "arm2" ]; then ARM2_GOLDEN="$2"; else ARM1_GOLDEN="$2"; fi
            shift 2 ;;
        --timeout) GEN_TIMEOUT="$2"; shift 2 ;;
        --hard-timeout) HARD_TIMEOUT="$2"; shift 2 ;;
        --shard) SHARD="$2"; shift 2 ;;
        --band) BAND="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        *) echo "error: unrecognized option $1" >&2; exit 1 ;;
    esac
done

if [ -z "$HARD_TIMEOUT" ]; then
    HARD_TIMEOUT="$(awk -v t="$GEN_TIMEOUT" 'BEGIN{printf "%d", t*1.5}')"
fi

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
    # 61 -> 62 in v0.4.10 with the golden re-freeze; see gate_arm1_encode.py's
    # EXPECTED_FIXTURE_COUNT for why the two counts are asserted independently.
    if [ "$n_done" -ne 62 ]; then
        echo "[gate/arm1] FAIL: #DONE $n_done, expected 62 -- short run, refusing to trust it" >&2
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

    # --band: restrict to one runtime stratum, using the golden's band column (7). Lets the
    # fast strata be gated routinely instead of only at release time -- the full v0.4.9
    # cohort is ~10 CPU-h serially, which nobody runs twice for a reproducibility check.
    if [ -n "$BAND" ]; then
        if [ ! -f "$ARM2_GOLDEN" ]; then
            echo "[gate/arm2] FAIL: --band needs a golden with a band column ($ARM2_GOLDEN missing)" >&2
            return 1
        fi
        names="$(awk -F'\t' -v b="$BAND" '$0 !~ /^#/ && $7 == b {print $1}' "$ARM2_GOLDEN" \
                 | grep -Fx -f <(echo "$names") || true)"
        echo "[gate/arm2] --band $BAND -> $(echo "$names" | grep -c .) molecule(s)" >&2
    fi

    # --shard i:n, ONE-BASED to match tools/test_dataset_roundtrip.py. `--shard 0:6` there
    # exits 2 with "index out of range", and someone launching 0..5 silently drops 1/6 of the
    # cohort. Same convention here so the two cannot be confused.
    local shard_i="" shard_n=""
    if [ -n "$SHARD" ]; then
        shard_i="${SHARD%%:*}"; shard_n="${SHARD##*:}"
        if ! [ "$shard_i" -ge 1 ] 2>/dev/null || [ "$shard_i" -gt "$shard_n" ]; then
            echo "error: --shard index out of range: '$SHARD' (1-based, expected 1:$shard_n .. $shard_n:$shard_n)" >&2
            return 2
        fi
        names="$(echo "$names" | awk -v i="$shard_i" -v n="$shard_n" 'NF && (NR-1)%n == i-1')"
        echo "[gate/arm2] shard $shard_i/$shard_n" >&2
    fi

    local expected
    expected="$(echo "$names" | grep -c .)"
    if [ "$expected" -eq 0 ]; then
        echo "[gate/arm2] FAIL: 0 molecules selected from $COHORT_DIR -- refusing an empty corpus" >&2
        return 1
    fi
    echo "[gate/arm2] cohort_dir=$COHORT_DIR molecules=$expected timeout=${GEN_TIMEOUT}s hard=${HARD_TIMEOUT}s out=$out" >&2
    echo "[gate/arm2] one subprocess per molecule -- this arm is EXPENSIVE by construction" >&2

    local i=0
    while IFS= read -r name; do
        [ -z "$name" ] && continue
        i=$((i + 1))
        echo "[gate/arm2] ($i/$expected) $name ..." >&2
        if [ "$HARD_TIMEOUT" -gt 0 ] 2>/dev/null; then
            timeout -s KILL "$HARD_TIMEOUT" "$PY" "$REPO/tools/gate_arm2_roundtrip_one.py" \
                --cohort-dir "$COHORT_DIR" --molecule "$name" --timeout "$GEN_TIMEOUT" \
                >> "$out" || true
        else
            "$PY" "$REPO/tools/gate_arm2_roundtrip_one.py" \
                --cohort-dir "$COHORT_DIR" --molecule "$name" --timeout "$GEN_TIMEOUT" \
                >> "$out"
        fi
        # A SIGKILLed child writes nothing, which would make #DONE short and fail the whole
        # arm as "truncated". Synthesize the row instead -- a hard timeout is a RESULT, and
        # the sentinel discipline this script documents says a missing line must never be
        # able to look like agreement. Guarded on absence so a child that printed and was
        # then killed while exiting cannot produce a duplicate.
        if ! grep -qE "^${name}"$'\t' "$out"; then
            printf '%s\tHARD_TIMEOUT@%ss\tHARD_TIMEOUT@%ss\t-\t-\t-\t\tERROR:HardTimeout:SIGKILL after %ss\n' \
                "$name" "$HARD_TIMEOUT" "$HARD_TIMEOUT" "$HARD_TIMEOUT" >> "$out"
        fi
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
    # (column 7 of the FRESH row, an observation-only column -- see
    # gate_arm2_roundtrip_one.py). Columns 7-8 of a v0.4.9 GOLDEN are band/honest_class,
    # also observation-only; nothing past column 3 is ever compared.
    #
    # SENTINELS (v0.4.9). A golden hash of the form NO_*@<budget>s is NOT a hash and must
    # never be compared as one:
    #   NO_STRUCTURE@Ns  the source run produced no smiles_2 at that budget. Whether it
    #                    produces one now is a fact about the BOX, not the code -- ULODUU
    #                    assembles at a 60 s cap and not at 30 s, which is precisely why the
    #                    boron fast-fail was refuted. Gate on sha1 only; report the outcome.
    #   NO_ENCODE@Ns     the ENCODER was hard-killed at the budget. No string exists at all,
    #                    so this row carries no byte-identity signal in either column.
    #   HARD_TIMEOUT@Ns  written by this script when it SIGKILLs a molecule (fresh side only).
    # Counting them separately is the point: "271 compared, 50 sha1-only, 4 observation" is
    # an honest verdict; folding them into a pass rate is how a gate becomes a clock.
    local mismatches=0 missing=0 compared=0 sha1_only=0 observed=0 det=0 not_run=0
    while IFS=$'\t' read -r g_name g_sha1 g_sha2 _rest; do
        [ -z "$g_name" ] && continue
        # A golden row for a molecule this invocation did not run -- because of --shard,
        # --band, or a cohort dir that is a subset -- is NOT missing. Only a molecule that
        # WAS in the run list and produced no row is missing, and that is a real failure.
        if ! echo "$names" | grep -qFx "$g_name"; then
            not_run=$((not_run + 1))
            continue
        fi
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

        case "$g_sha1" in
            NO_ENCODE@*)
                observed=$((observed + 1))
                echo "[gate/arm2] OBSERVE $g_name (golden $g_sha1): fresh sha1=${f_sha1:-<none>} sha2=${f_sha2:-<none>}" >&2
                continue ;;
        esac

        if [ "$f_sha1" != "$g_sha1" ]; then
            echo "[gate/arm2] MISMATCH: $g_name (smiles_1)" >&2
            echo "[gate/arm2]   golden sha1=$g_sha1" >&2
            echo "[gate/arm2]   fresh  sha1=$f_sha1" >&2
            echo "[gate/arm2]   fresh row: $row" >&2
            mismatches=$((mismatches + 1))
            continue
        fi

        case "$g_sha2" in
            NO_STRUCTURE@*)
                sha1_only=$((sha1_only + 1))
                if [ -n "$f_sha2" ]; then
                    echo "[gate/arm2] OBSERVE $g_name: produced no structure at the source budget, produces one NOW (sha2=$f_sha2)" >&2
                else
                    echo "[gate/arm2] OBSERVE $g_name: still produces no structure" >&2
                fi
                continue ;;
            NO_STRUCTURE_DET)
                # A deterministic failure IS a code property -- gate it.
                if [ -n "$f_sha2" ]; then
                    echo "[gate/arm2] MISMATCH: $g_name deterministically produced no structure in the" >&2
                    echo "[gate/arm2]   source run and produces one now (sha2=$f_sha2). If that is intended," >&2
                    echo "[gate/arm2]   rebuild the golden and say so." >&2
                    mismatches=$((mismatches + 1))
                    continue
                fi
                det=$((det + 1))
                continue ;;
        esac

        if [ "$f_sha2" != "$g_sha2" ]; then
            echo "[gate/arm2] MISMATCH: $g_name (smiles_2)" >&2
            echo "[gate/arm2]   golden sha2=$g_sha2" >&2
            echo "[gate/arm2]   fresh  sha2=$f_sha2" >&2
            echo "[gate/arm2]   fresh row: $row" >&2
            mismatches=$((mismatches + 1))
            continue
        fi
        compared=$((compared + 1))
    done < <(grep -vE '^#' "$ARM2_GOLDEN")

    echo "[gate/arm2] verdict basis: $compared fully compared, $det deterministic-no-structure," >&2
    echo "[gate/arm2]   $sha1_only sha1-only (NO_STRUCTURE, budget-dependent), $observed observation-only (NO_ENCODE)," >&2
    echo "[gate/arm2]   $not_run golden row(s) not in this run's selection" >&2
    if [ $mismatches -gt 0 ] || [ $missing -gt 0 ]; then
        echo "[gate/arm2] FAIL: $mismatches mismatch(es), $missing missing" >&2
        return 1
    fi
    if [ $((compared + det)) -eq 0 ]; then
        echo "[gate/arm2] FAIL: 0 molecules were actually gated -- a green verdict on an" >&2
        echo "[gate/arm2]       all-sentinel selection certifies nothing. Refusing to pass." >&2
        return 1
    fi
    echo "[gate/arm2] PASS -- $((compared + det)) gated (+$sha1_only sha1-only, $observed observed)" >&2
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
        echo "usage: $0 <arm1|arm2|both> [--fixtures-dir DIR] [--cohort-dir DIR] [--golden PATH]" >&2
        echo "         [--timeout S] [--hard-timeout S] [--shard i:n] [--band NAME] [--out PATH]" >&2
        echo "  --shard is ONE-BASED (1:6 .. 6:6), matching tools/test_dataset_roundtrip.py" >&2
        echo "  --band  filters arm2 to one runtime stratum from the golden's band column" >&2
        exit 1
        ;;
esac
