# Process Document: Bidentate Placement Fidelity (TASK-30/31/32)

**Generated:** 2026-07-03T18:01:45-07:00
**Session Focus:** Diagnose and fix the placement-geometry bug that broke the
byte-stable XYZ→OIN→XYZ round trip for bidentate chelates (DIPAMP), and clear
the regression that fix surfaced in the Zone-A P stereo-enforcement tests.

## Problem Statement

Two round-trip tests were `@unittest.expectedFailure`:
`test_p_stereocenter_roundtrip` (DIPAMP, a bidentate P^P chelate on Rh) and
`test_haptic_face_golden_match` (an eta-ring case, out of scope here). A prior
session had already pinned the DIPAMP root cause to a placement bug in
`_stitch_fragment` (`src/oinsmiles/generation/molassembler_adapter.py`): the
bidentate Kabsch-alignment path has one rotational degree of freedom (around
the bite axis) that a bite-axis optimizer resolves, but the optimizer's safety
guard only checks non-binding **heavy** atoms for metal collisions, not
hydrogens. The goal of this arc was to measure precisely which failure
mechanism was at play — a fixable objective/guard bug, or a deeper conformation
mismatch requiring a re-embed redesign — and fix it without guessing.

## Starting State

Repo at commit `aefd75f` (worklog entry for TASK-30 already added, no `src/`
changes yet). Full unit suite green except two `expectedFailure` xfails.
`src/oinsmiles/generation/molassembler_adapter.py::_stitch_fragment` had:
- a bite-distance compatibility check (`:1493-1498`) rejecting only when the
  isolated ligand's binding-atom distance differs from the chelate target by
  more than 2.0 Å;
- a 360°/5°-step bite-axis rotation optimizer (`:1511-1571`) that already
  factors H atoms into its own minimize-collision objective;
- a post-placement rejection guard (`:1573-1584`) that only inspects
  non-binding **heavy** atoms (`symbols[i] != "H"`) against a 1.7 Å threshold.

`tests/unit/test_zone_a_p_genenforce.py` (Phase 4b, done earlier the same
session) used the DIPAMP fixture to exercise Zone-A P stereo enforcement on
the Kabsch/template placement path, via a `_test_flip_chiral_idx` forced
mis-embed hook.

## Approach & Methodology

Diagnostic-first, matching the pattern established earlier in this session
(TASK-04, TASK-10, TASK-20): write a throwaway measurement script, decide the
fix class from evidence, *then* write the smallest fix that matches the
evidence, then re-run the full suite to catch anything the fix disturbs. Three
lightweight (non-HACF) tasks in sequence — TASK-30 (diagnose) → TASK-31 (fix)
→ TASK-32 (repair a fix-induced regression) — rather than a HACF MiniPRD,
because the fix stayed small enough that the task-file pattern already in use
this session was sufficient.

## Steps Taken

1. **TASK-30 — instrumented diagnostic (Sonnet).** Wrote a scratch script
   (job scratchpad, not committed) that: (a) ran the real, unmodified
   XYZ→OIN→XYZ pipeline on `Rh-RR-DIPAMP-Cl2.xyz` and `fac_irppy3.xyz` and
   reported per-element nearest-neighbour-to-metal distances; (b) monkeypatched
   `molassembler_adapter._stitch_fragment` to *capture* the real call's
   arguments (fragment SMILES, binding indices, target positions, forbidden
   positions, seed) without altering its behavior, then replayed the bidentate
   branch offline with the same ETKDG seed (deterministic conformer) to log
   the **full** 72-angle sweep instead of only the one angle production code
   keeps; (c) forced the guard to reject DIPAMP's placement (in-process
   monkeypatch only) to see what the existing DG fallback produces.
   Result: DIPAMP's isolated ligand P···P = 4.408 Å vs. chelated target
   3.182 Å (1.226 Å mismatch, inside the 2.0 Å compatibility threshold so it's
   accepted into the Kabsch path); the production algorithm's chosen angle
   (100°) already equals the sweep's **global optimum** for H-to-metal
   distance (1.390 Å) — 0 of 72 angles reach a safe distance, worst case
   0.261 Å. This ruled out an objective/grid bug: no rotation fixes it. Forcing
   DG fallback produced a clean chelate (H min 3.19 Å, Rh–P 2.41/2.42 Å) that
   re-encoded to a **byte-identical** OIN. Decision recorded: root cause is a
   genuine conformation mismatch (mechanism "B"), but the fix is cheap because
   DG already handles it correctly — the only blocker is the heavy-atom-only
   guard. No `src/` changed this step; `spec/worklog/TASK-30-*.md` `Status:`
   and `NOTES.md` updated with the full numbers.

2. **TASK-31 — implement the guard fix (Sonnet).** Added a second check
   immediately next to the existing `:1573-1584` heavy-atom guard: for
   bidentate+ placements (`len(binding_idxs) >= 2`), also reject (return
   `None` → DG fallback) if any non-binding **H** atom lands within 1.8 Å of
   the metal centre. Left the heavy-atom check, the bite-axis optimizer, and
   the DG worker untouched, per the task's explicit constraints (do not
   tighten the `:1497` bite-delta threshold as a proxy — an outcome-based H
   check was preferred so clean bidentates with unknown deltas aren't put at
   risk). Verified: the real pipeline now round-trips DIPAMP byte-stable, and
   flipped `test_p_stereocenter_roundtrip`'s `expectedFailure` to a hard pass.
   Running the full suite surfaced an **unanticipated regression**: 2 failures
   + 1 error in `tests/unit/test_zone_a_p_genenforce.py` (Phase 4b). Root
   cause traced: those 3 tests depended on DIPAMP's P^P fragment reaching the
   Kabsch/template path (to exercise the mis-embed-detection/re-embed-retry
   machinery via `_test_flip_chiral_idx`) — but the new guard now
   *unconditionally* routes DIPAMP to DG (the 1.226 Å bite mismatch is a fixed
   geometric fact of this fixture, present on every attempt regardless of
   chirality flips), so that machinery becomes permanently unreachable for
   this fixture. Correctly identified as **out of scope** for TASK-31 (whose
   constraints named only `molassembler_adapter.py` and one test file) and
   left `Status: BLOCKED` rather than force a cross-file fix — changes kept
   uncommitted for review.

3. **TASK-32 — retarget the blocked tests (Sonnet).** Rather than touch the
   guard or the placement logic again, retargeted the 3 affected
   `test_zone_a_p_genenforce.py` tests off DIPAMP onto synthetic monodentate
   P-stereocenter OIN strings that verifiably stay on the template/Kabsch
   path (`_MONO_P_OIN` on `[Ni_TET]`, plus `_MONO_P_CORESIDENT_OIN` with a
   directly-bonded carbon stereocenter for the non-interference check) —
   chosen because Zone-A P enforcement (`_verify_zone_a_p`) is denticity-
   independent; a monodentate P exercises the identical mechanism. Both
   fixtures were confirmed empirically (scratch one-liners) *before* being
   wired into the test file, per the session's established discipline of
   verifying before committing to a test edit. Added a new test,
   `test_dipamp_dg_fallback_warns_honestly`, documenting TASK-31's routing
   decision as an intentional, warned trade-off (DIPAMP legitimately takes
   the DG path now and must emit `OINStereoWarning`, not silently pass).
   Result: `test_zone_a_p_genenforce.py` back to 12/12 passing; full unit
   suite 113 run, OK (skipped=3, expected failures=1 — only
   `test_haptic_face_golden_match` remains); root suite 55 run, OK.
   TASK-31 flipped `BLOCKED` → `DONE` (blocker cleared); TASK-32 → `DONE`.

4. **Follow-on discovery (same session, informs but is out of scope for this
   arc):** while verifying no other regressions, an enforcement asymmetry
   was noticed and confirmed — `[P@]` on square-planar (`SPL`) complexes
   enforces to the wrong enantiomer (the re-embed-only correction mechanism
   can't reach it, since the metal-present CIP on SPL is fixed by which
   pyramid face the metal lands on, not by conformer resampling). This was
   spun out as its own workstream (`SPL-P-enforcement-design-brief.md` →
   `/hyper-consult-cto` → `SPL-P-enforcement-decision.md` → MiniPRD-C via
   `/hyper-architect`/`/hyper-redteam`/`/hyper-resolve`) and is documented in
   the separate process docs `process_20260703_171405_spl-zone-a-p-
   enforcement-decision.md` and `process_20260703_173828_zone-a-p-spl-
   dummy-embed-resolve.md` — not repeated here.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| Root-cause mechanism = conformation mismatch ("B"), not objective/grid bug ("A") | Assuming "A" and just widening the bite-axis search grid or reweighting its objective | Empirically ruled out: the production algorithm's chosen angle already equals the full-360°-sweep global optimum for H-to-metal distance; 0/72 angles are safe |
| Fix = outcome-based non-binding-H guard (reject if H < 1.8 Å from metal) | Tightening the `:1497` bite-delta compatibility threshold instead | The delta is a proxy with unknown safe range for other, currently-clean bidentates; the H-distance outcome is what actually causes the round-trip corruption, so gating on it directly is lower-risk |
| Retarget the 3 blocked Phase-4b tests to a synthetic monodentate fixture rather than relax TASK-31's guard | Loosening the new H-guard so DIPAMP could still reach the Kabsch path some of the time | Zone-A P enforcement is denticity-independent, so a monodentate fixture exercises the identical mechanism without reopening the placement bug TASK-31 just closed |
| Stop and mark TASK-31 `BLOCKED` rather than fix the regression inline | Extending TASK-31's own scope to also patch `test_zone_a_p_genenforce.py` | The regression lived in a different file than TASK-31's stated constraints permitted; treated as a scope boundary requiring a follow-on task (TASK-32), not a unilateral expansion |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Bite-axis diagnostic script | job scratchpad (not in repo) | created, discarded after use |
| TASK-30 spec | `spec/worklog/TASK-30-bidentate-placement-diagnostic.md` | `Status:` updated with decision summary |
| TASK-31 spec | `spec/worklog/TASK-31-bidentate-placement-fix.md` | created, then `Status: DONE` after TASK-32 cleared its blocker |
| TASK-32 spec | `spec/worklog/TASK-32-retarget-phase4b-kabsch-tests.md` | created, `Status: DONE` |
| Bidentate guard fix | `src/oinsmiles/generation/molassembler_adapter.py` | added non-binding-H rejection check next to the existing heavy-atom guard (`_stitch_fragment`) |
| DIPAMP round-trip test | `tests/unit/test_stereo_roundtrip_diagnostics.py` | removed `expectedFailure` from `test_p_stereocenter_roundtrip`; docstring updated |
| Zone-A P enforcement tests | `tests/unit/test_zone_a_p_genenforce.py` | 3 tests retargeted from DIPAMP to synthetic monodentate fixtures; 1 new test (`test_dipamp_dg_fallback_warns_honestly`) added |
| Session worklog | `spec/worklog/NOTES.md` | dated Log entries for TASK-30, TASK-31 (incl. the blocked regression), TASK-32; status table rows updated |

## Results & Outcomes

- `test_p_stereocenter_roundtrip` is now a hard pass — DIPAMP round-trips
  XYZ→OIN→XYZ→OIN byte-identically, including its `[P@]{0}`/`[P@]{1}` tags,
  via the DG fallback path.
- `uv run python -m unittest discover tests/unit` → 113 run, OK (skipped=3,
  expected failures=1 — only `test_haptic_face_golden_match`, an unrelated
  eta-ring case, remains xfail).
- `uv run python -m unittest discover tests` → 55 run, OK.
- No regression on previously-clean generators: cisplatin, transplatin,
  cis-PtCl2(en), fac/mer-Ir(ppy)3 all still produce their expected strings
  (confirmed via `tests.unit.test_regression_stability`, re-run after TASK-32).
- Uncovered and separately scoped a real correctness bug (SPL Zone-A P
  enforcement is one-sided) as a follow-on workstream, now past its
  `/hyper-consult-cto` decision and `/hyper-resolve`d as MiniPRD-C, ready for
  `/hyper-execute` (tracked in `spec/worklog/NOTES.md` and its own process
  docs, not this one).
- As of this writing, the TASK-31 guard fix and TASK-32 test retargets remain
  **uncommitted**, left staged for human review per the tasks' own
  instructions (this session does not commit source-affecting changes
  without being asked).

## How to Reproduce

Starting from a clean checkout at commit `aefd75f` (or later, before `ee0b3f0`):

1. Read `spec/worklog/TASK-30-bidentate-placement-diagnostic.md` and follow it
   exactly: build a throwaway diagnostic script (do not touch `src/`) that
   captures `_stitch_fragment`'s real call arguments via monkeypatch, replays
   its bidentate bite-axis sweep with full 72-angle instrumentation, and
   separately forces the DG-fallback path to check its output quality. Record
   the measured numbers and the A/B/mix decision in `spec/worklog/NOTES.md`
   and the task file's `Status:` line.
2. Read `spec/worklog/TASK-31-bidentate-placement-fix.md` and implement
   exactly the described guard addition in
   `src/oinsmiles/generation/molassembler_adapter.py::_stitch_fragment`
   (non-binding H within 1.8 Å of the metal → reject, bidentate+ only), then
   flip `test_p_stereocenter_roundtrip`'s `expectedFailure`. Run
   `uv run python -m unittest discover tests/unit` and `discover tests`. If a
   regression appears outside the two files TASK-31 names, stop and mark the
   task `BLOCKED` rather than expanding scope — write up the exact regression
   and root cause in `NOTES.md`.
3. Read `spec/worklog/TASK-32-retarget-phase4b-kabsch-tests.md` (or design an
   equivalent fixture-retarget) to move any tests that depended on the
   now-DG-routed fixture onto a synthetic fixture that verifiably stays on the
   affected code path. Verify each replacement fixture empirically (a scratch
   check) before wiring it into the test file. Re-run the full suite and
   confirm both blocked tasks clear.

Gotcha: the guard fix and the Phase-4b test dependency are coupled through the
*same underlying fixture* (DIPAMP) being used for two unrelated things (a
placement-geometry regression test and a stereo-enforcement mechanism test) —
fixing one without checking the other's assumptions is what caused the
TASK-31 blocker. Any future fix that changes which code path a fixture
resolves to should grep test files for other uses of that fixture first.

## Patterns & Lessons

- **Measure before fixing, even when the fix looks obvious.** The bite-axis
  guard's "only checks heavy atoms" gap looked like a simple oversight, but
  the diagnostic's real value was ruling out the more expensive hypothesis
  (an objective/grid-search bug needing a redesigned optimizer) — without
  that step, TASK-31 might have chased a rotation-search fix that could never
  have worked, since no angle is actually safe for this fixture.
  See [v3.7 worklog process](../worklog/NOTES.md) for the diagnostic-first
  convention this session established and reused across TASK-04/10/20/30.
- **A shared test fixture across two unrelated test suites is a hidden
  coupling.** `Rh-RR-DIPAMP-Cl2.xyz` was independently useful for (a) a
  placement-geometry round-trip test and (b) a stereo-enforcement mechanism
  test — a fix that changes fixture (a)'s code path silently broke (b).
  When retargeting or fixing behavior tied to a specific fixture, grep the
  whole test tree for other consumers of that fixture before declaring done.
- **Stopping at a scope boundary and filing a follow-on task (rather than
  quietly widening the current task) kept the acceptance trail honest** —
  TASK-31's `BLOCKED` status accurately reflected that 3 of its 4 acceptance
  checks passed and one didn't, rather than silently absorbing an unrelated
  file into its diff.
