#!/bin/bash
# v0.4.5 live state — run this FIRST in a fresh session, before reading anything else.
#
# Prints branch tips, uncommitted work, verification status and the next action, computed
# from the repo rather than from a doc, so it cannot go stale. Written after an account
# spend limit killed eight parallel lane agents mid-task; the lesson was that a
# hand-maintained status file is stale the moment work resumes.
#
#   tools/v045_state.sh            # full state
#   tools/v045_state.sh --next     # just the next action
set -u
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1

LANES="1 2 3 4 5 6 7 8 perf"
BAR="────────────────────────────────────────────────────────────────────────"

branch_for() { [ "$1" = perf ] && echo "swimlane/v045-perf" || echo "swimlane/v045-lane$1"; }
tree_for()   { [ "$1" = perf ] && echo "../oin-v045-perf"  || echo "../oin-v045-lane$1"; }

if [ "${1:-}" != "--next" ]; then
echo "$BAR"
echo " v0.4.5 STATE  @ $(date -Is)"
echo "$BAR"
echo
echo "MAIN (local, unpushed — standing instruction: do NOT push)"
echo "  tip: $(git log --oneline -1 --format='%h %s' main)"
echo "  commits since the v0.4.5 base e8b603d5: $(git rev-list --count e8b603d5..main)"
echo "  encoder defaults changed on main: NONE (every new lever is default-OFF)"
echo
echo "LANES"
printf "  %-6s %-22s %6s %5s  %s\n" lane branch ahead WIP tip
for L in $LANES; do
    BR=$(branch_for "$L"); WT=$(tree_for "$L")
    if git rev-parse --verify -q "$BR" >/dev/null 2>&1; then
        AHEAD=$(git rev-list --count "main..$BR" 2>/dev/null || echo "?")
        DIRTY=$(git -C "$WT" status --porcelain 2>/dev/null | wc -l)
        TIP=$(git log --oneline -1 --format='%s' "$BR" | cut -c1-40)
        MARK=""
        git log -1 --format='%s' "$BR" | grep -q '^WIP' && MARK=" <-- WIP/INCOMPLETE"
        printf "  %-6s %-22s %6s %5s  %s%s\n" "$L" "${BR#swimlane/}" "$AHEAD" "$DIRTY" "$TIP" "$MARK"
    else
        printf "  %-6s %-22s %6s %5s  %s\n" "$L" "-" "-" "-" "NEVER STARTED"
    fi
done
echo
echo "  Uncommitted files in a lane worktree are FRAGILE. Commit them as"
echo "  'WIP(laneN): ... INCOMPLETE' before doing anything else."
echo
echo "LEVERS introduced by v0.4.5 (all default OFF; nothing promoted yet)"
echo "  searched across lane branches, since most are not on main yet:"
for L in $LANES; do
    BR=$(branch_for "$L")
    git rev-parse --verify -q "$BR" >/dev/null 2>&1 || continue
    for f in $(git ls-tree -r --name-only "$BR" -- src/ | grep '\.py$'); do
        git show "$BR:$f" 2>/dev/null \
          | grep -ohE 'OIN_(CANONICAL|STABLE)_[A-Z_]+' || true
    done | sort -u | sed "s/^/    lane$L: /"
done
echo "    (plus OIN_EMIT_AXIAL, pre-existing, still default OFF per the product call)"
echo
echo "LONG JOBS (survive a harness timeout; poll instead of holding a task open)"
systemctl --user list-units --state=active 'lane*' 'v045*' --no-legend 2>/dev/null \
  | awk '{print "  "$1" "$3"/"$4}' || true
[ -z "$(systemctl --user list-units --state=active 'lane*' 'v045*' --no-legend 2>/dev/null)" ] \
  && echo "  (none active)"
echo "  load: $(uptime | sed 's/.*load average: //')   cores: $(nproc)"
echo
echo "READ THESE, IN ORDER"
for d in docs/V045_STATUS_2026-07-25.md docs/RENUMBERING_INSTABILITY_v0.4.5.md \
         docs/CANONICAL_BODY_v0.4.5.md docs/WINDING_RESIDUAL_v0.4.5.md \
         docs/AXIAL_v0.4.5_LANE4.md; do
    [ -f "$d" ] && echo "  $d" || echo "  $d   (on a lane branch, not main)"
done
echo
fi

echo "NEXT ACTION"
NEXT=""
for L in 2 8; do
    BR=$(branch_for "$L")
    git rev-parse --verify -q "$BR" >/dev/null 2>&1 || continue
    if git log -1 --format='%s' "$BR" | grep -q '^WIP'; then
        [ -z "$NEXT" ] && NEXT="$L"
    fi
done
if [ "$NEXT" = 2 ]; then
    cat <<'EOT'
  Finish Lane 2 (CRITICAL PATH — Lanes 5 and 6 are blocked on it).
    worktree ../oin-v045-lane2, branch swimlane/v045-lane2
    1. DUDREA_comp_0 is DIAGNOSED and FIXED behind OIN_STABLE_METAL_AC (commit
       8bf9df61). It was never a slot-labelling problem: xyz2AC_obabel's
       valence-capping loop iterated in atom-index order, so capping the metal
       before vs after a bridging hydride decided whether the Y-H bond survived
       (degree 5 / SPY vs degree 4 / TET). AC differed in 3/8 renumberings off,
       0/8 on. STILL NEEDED: a corpus A/B, because this changes PERCEPTION and
       capping the metal first can only ADD bonds the old order discarded.
    2. Lane 2's own slot post-pass still has NO A/B numbers. Measure with:
         PYTHONPATH=$PWD/src <venv>/python tools/canonicality_probe.py --n 200 --trials 2
       once with levers off, once with OIN_CANONICAL_SLOTS=1.
    3. Over-folding guards NOT confirmed green — these are the whole risk:
         tests/unit/test_facmer_key.py  tests/integration/test_isomer_divergence.py
EOT
elif [ "$NEXT" = 8 ]; then
    echo "  Finish Lane 8: extend oin/stable_stereo.py to trivalent P donors, then"
    echo "  prove the descriptor still FLIPS for the mirror (not stable-because-constant)."
else
    cat <<'EOT'
  INTEGRATION IS DONE. All 16 lanes are merged on `release/v0.4.5` (worktree
  ../oin-v045-final), 84 commits ahead of main, ruff clean, and the six canonicality
  levers are PROMOTED to default-ON through src/oinsmiles/oin/levers.py.

  Remaining, in order:
    1. Triage the release suite. Promotion makes every test written as "levers OFF ->
       byte-identical" fail BY DESIGN -- Lane 2 established that CisPlatin, TransPlatin
       and fac-Ir(ppy)3 move under OIN_CANONICAL_SLOTS (key-identical, pinned separately
       in TestLeverOnGoldens). Separate "asserts the old default" from real breakage and
       update the former. Do NOT merge to main until this is clean.
    2. Read the final re-baseline number once tier 2 drains. Check
       `systemctl --user is-active v045-rebaseline` -- NOT the report count. 936/936
       reports existed while the g-xTB tier was still running, and reading count as
       completion once already produced a premature published figure.
    3. Merge release/v0.4.5 -> main (local only, DO NOT PUSH), tag before deleting any
       swimlane branch.
    4. Optional: the 5k sweep on cohort-v0.4.5-5k, now that the config is settled.

  Two levers deliberately held OFF with reasons in levers.py::_HELD_OFF -- read those
  before promoting either: OIN_H_FAITHFUL is BLOCKED (canonical_body_emit's reparse undoes
  it) and OIN_BORON_CAGE is decided-to-promote but not yet enabled (it moves 14
  scored-passing molecules to failing, which is correct -- they were wrong).
EOT
fi
echo
echo "RULES THAT BIT US ALREADY"
cat <<'EOT'
  * Use the MAIN checkout's .venv (rdkit pinned 2025.9.3). NEVER uv sync in a worktree.
  * Tests: `discover tests/unit`, NOT `discover tests`.
  * Commit normally, NEVER --no-verify; include a Claude-Session: trailer (the hook
    rewrites it; omit it and you get no trailer at all).
  * NEVER `git stash` — it is shared across worktrees and a sibling will collide.
    A/B by `git show HEAD:path > path` with an EXIT trap instead.
  * Wall-clock timing is meaningless above ~load 12. Measure deterministic counters.
  * Long jobs: systemd-run --user -p OOMPolicy=continue -p MemoryMax=8G, then poll
    `systemctl --user is-active`. A harness task timeout WILL kill a plain background job.
  * Prefer inline work over spawning subagents: eight parallel agents at ~350k tokens
    each is what hit the spend limit.
EOT
