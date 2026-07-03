# TASK-31: Route incompatible-bite chelates to DG fallback (placement fix)

Status: DONE — guard change + xfail flip both implemented and verified
correct (acceptance #1 and #2 both PASS: real-pipeline DIPAMP round-trip is
byte-stable, `test_p_stereocenter_roundtrip` passes). The blocker on
acceptance #3 is CLEARED by TASK-32: `tests/unit/test_zone_a_p_genenforce.py`
(Phase 4b GenEnforce)'s 3 tests that depended on DIPAMP's P^P fragment
reaching the Kabsch/template placement path have been retargeted to a
synthetic monodentate P-stereocenter fixture that verifiably stays on that
path (see `spec/worklog/TASK-32-retarget-phase4b-kabsch-tests.md`, DONE).
Full suite now green end-to-end: unit skipped=3/expected failures=1
(only `test_haptic_face_golden_match`), root suite OK,
`test_regression_stability` OK. See the dated 2026-07-03 Log entries in
`spec/worklog/NOTES.md` (TASK-31 and TASK-32) for full acceptance output and
root-cause detail. Not committed; left unstaged for review per instructions.
Depends on: TASK-30 (diagnostic — DONE, mechanism = B/conformation, DG fallback
already produces clean output)
Suggested model: Sonnet

## Goal

Make DIPAMP-class bidentate ligands (whose isolated conformation can't span the
chelate bite) fall back to the DG builder, which already generates them
correctly. TASK-30 proved: for `Rh-RR-DIPAMP-Cl2.xyz` no rigid bite-axis
rotation avoids the metal collision (0/72 angles safe; chosen angle is already
the global optimum), but forcing the existing DG fallback yields a clean
chelate (Rh–P 2.41 Å, min non-binding H 3.19 Å) AND a **byte-identical**
round-trip. The only thing blocking that today is a guard that misses the
collision.

## Root cause (exact, from TASK-30)

`src/oinsmiles/generation/molassembler_adapter.py`, `_stitch_fragment`:
- `:1493-1498` bite-distance compatibility check rejects only when
  `abs(org_dist - tgt_dist) > 2.0`. DIPAMP's delta is **1.226 Å** (isolated
  P···P 4.408 Å vs target 3.182 Å) → passes → enters Kabsch path.
- `:1573-1584` post-placement rejection guard checks only non-binding **heavy**
  atoms (`symbols[i] != "H"`) against 1.7 Å. DIPAMP's worst heavy atom is
  1.754 Å (clears by 0.054 Å), while 3 non-binding **H** atoms land at
  1.39–1.65 Å from Rh — undetected. Those H atoms are then perceived as Rh–H
  hydrides by XYZ→OIN, corrupting the round-trip topology.

## The fix (primary — regression-safe by construction)

Extend the `:1573-1584` guard to also reject on non-binding **H** atoms near
the metal. Rationale for safety: this branch is bidentate/polydentate only
(`len(binding_idxs) >= 2`) and checks NON-binding atoms; a genuine metal
hydride's H is a BINDING atom, so it is excluded — hydride complexes are
unaffected. A non-binding H within ~1.8 Å of a metal is always a spurious
clash (clean chelates keep non-binding H ≳ 2.5 Å; DIPAMP-via-DG gives 3.19 Å).

Concretely: add a non-binding-H check (suggest threshold **1.8 Å**; keep the
existing heavy-atom 1.7 Å check) that returns `None` → DG fallback. Keep it a
small, readable addition next to the existing guard; do not rewrite the guard.

Do NOT instead just tighten the `:1497` bite-delta threshold — that is a proxy
(we don't know clean bidentates' deltas, so it risks regressing them). The
outcome-based H guard is preferred. (You may additionally note the delta in a
comment, but the H guard is the mechanism.)

## Acceptance (exact commands + expected)

1. Real end-to-end round-trip is now byte-stable (independently re-verify
   TASK-30's monkeypatched claim on the REAL pipeline):
   ```
   uv run python -c "
   from oinsmiles import XYZToSMILES
   from oinsmiles.generation.engine import OIN3DGenerator
   import tempfile, os, warnings
   warnings.simplefilter('ignore')
   o1 = XYZToSMILES().convert('tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz')
   g = OIN3DGenerator().generate(o1)
   with tempfile.NamedTemporaryFile('w', suffix='.xyz', delete=False) as fh:
       fh.write(g.xyz); p = fh.name
   o2 = XYZToSMILES().convert(p); os.unlink(p)
   assert o1 == o2, f'NOT byte-stable:\n{o1}\n{o2}'
   print('DIPAMP round-trip byte-stable:', o1)
   "
   ```
   Expect: prints the byte-stable OIN (with `[P@]{0}`/`[P@]{1}`).
2. Flip the now-passing xfail: in
   `tests/unit/test_stereo_roundtrip_diagnostics.py`, remove
   `@unittest.expectedFailure` from `test_p_stereocenter_roundtrip` (it should
   now pass). Update its docstring/comment to note it passes as of this fix via
   DG fallback. If it does NOT pass, STOP and report — do not force it.
3. Full suite green, no regression:
   ```
   uv run python -m unittest discover tests/unit 2>&1 | tail -3
   uv run python -m unittest discover tests 2>&1 | tail -3
   ```
   Expect unit: one fewer expected-failure than before (was skipped=3,
   expected failures=2 → now expected failures=1), everything else unchanged
   (`test_haptic_face_golden_match` stays xfail — it's an eta ring, a different
   path this fix does not touch).
4. No regression of currently-clean generators — confirm these still produce
   their expected strings (they must NOT newly route to DG):
   ```
   uv run python -m unittest tests.unit.test_regression_stability -v 2>&1 | tail -4
   ```
   Expect: cisplatin, transplatin, cis-PtCl2(en), fac/mer-Ir(ppy)3 all OK.

## Constraints / DO NOT

- Touch only `src/oinsmiles/generation/molassembler_adapter.py` (the guard) and
  `tests/unit/test_stereo_roundtrip_diagnostics.py` (un-xfail one test).
- Do NOT tighten the `:1497` bite-delta threshold as the mechanism.
- Do NOT alter the heavy-atom guard, the bite-axis optimizer, or the DG worker.
- Do NOT touch `test_haptic_face_golden_match` (separate eta-ring problem).
- Pre-commit hook is green now: run `uv run ruff format` on the file you edit
  before committing so the hook passes without `--no-verify`.

## On completion

Set `Status: DONE`, append a dated Log entry to `spec/worklog/NOTES.md`
(the guard change, the flipped xfail, acceptance results), update the TASK-31
row + the two-xfail note in the Session-state block (now one xfail:
`test_haptic_face_golden_match`). Do NOT commit — leave staged for review.
