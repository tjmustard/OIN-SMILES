# TASK-32: Retarget Phase-4b Kabsch-enforcement tests off DIPAMP

Status: DONE
Depends on: TASK-31 (placement fix — in the working tree, UNSTAGED; keep it)
Suggested model: Sonnet

## Why

TASK-31's placement fix routes DIPAMP to the DG fallback (its bite is
incompatible with the Kabsch path — proven in TASK-30). That is correct and
DIPAMP now round-trips losslessly. But 3 tests in
`tests/unit/test_zone_a_p_genenforce.py` used DIPAMP to exercise the
**template/Kabsch-path P-stereo enforcement machinery**, which DIPAMP no longer
touches. Retarget them to a P-stereocenter that STAYS on the template path.
Keep the TASK-31 fix; do not revert it.

## Verified facts (use these — measured 2026-07-03, do not re-derive)

- P-stereo enforcement runs at the ASSEMBLY level in `_template_generate`
  (`molassembler_adapter.py` `_verify_zone_a_p`, ~:1928/:1968), for ANY fragment
  carrying a Zone-A P tag that reaches the template path — denticity-independent.
  The forced-flip hook `_test_flip_chiral_idx` is per-fragment in
  `_stitch_fragment` (:1326/:1399). So a MONODENTATE P-stereocenter on the
  template path exercises the same machinery.
- Expected labels come from the OIN's own `[P@]`/`[P@@]` tag
  (`_zone_a_p_expected_labels`, :341) — so a SYNTHETIC hand-written OIN provides
  its own oracle. **No new XYZ fixture is required**; construct inline OINs like
  the existing fallback tests already do (`c1ccccc1[P@](CC)C`).
- **Verified clean-enforcing template-path fixtures** (mol built, 0
  OINStereoWarning): `[Ni_TET].c1ccccc1[P@]{0}(CC)C.[Cl]{1}.[Cl]{2}.[Cl]{3}`
  and `[Pt_TET].c1ccccc1[P@]{0}(CC)C.[Cl]{1}.[Cl]{2}.N{3}`. Use a TET geometry.
- **AVOID** `[P@]` on `[Pt_SPL]` — it warns "could not be enforced" (a noted
  enforcement asymmetry: `[P@@]` on the same SPL complex is clean). Pick a
  tag+geometry combo that enforces cleanly; verify each fixture you choose
  actually yields the warning count the test expects before wiring it in.

## Work

Retarget these 3 tests (in `tests/unit/test_zone_a_p_genenforce.py`); preserve
each test's ORIGINAL intent, only change the fixture/vehicle:

1. `test_dipamp_generation_is_clean` (TestZoneAPNoRegression) — intent: a
   P-stereocenter complex generates on the enforcing template path with 0
   OINStereoWarnings. Retarget to a verified clean-enforcing TET P fixture
   above. Rename to reflect it's a generic Kabsch-path P case (e.g.
   `test_monodentate_p_generation_is_clean`). SEPARATELY, DIPAMP now legitimately
   routes to DG and warns — if you want to keep a DIPAMP entry, assert it emits
   the honest DG "unenforced on fallback path" warning (do NOT assert silence).
2. `test_persistent_mismatch_warns_once_and_completes_quickly`
   (TestZoneAPBoundedFailure) — intent: a forced, un-correctable flip on the
   template path yields exactly ONE "could not be enforced" warning and still
   emits a mol, bounded in time. Retarget the vehicle to a template-path P
   fixture; keep the `_test_flip_chiral_idx` forced-flip mechanism and the
   single-"could not be enforced"-warning assertion.
3. `test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident`
   (TestZoneAPForcedMisEmbedCorrection) — intent: a forced single-atom
   mis-embed of a Zone-A P is corrected (0 warnings) WITHOUT disturbing a
   CO-RESIDENT stereocenter in the same fragment. DIPAMP supplied the co-resident
   via its 2nd P. Replace with a single ligand that has BOTH a Zone-A P
   stereocenter AND a co-resident stereocenter (e.g. a phosphine bearing a
   carbon stereocenter substituent, `...[P@](...)[C@H](C)...`). Construct it,
   confirm it stays on the template path and enforces cleanly at baseline, then
   apply the same forced-mis-embed + correction assertion. If no such co-resident
   fixture can be made to enforce cleanly, STOP and report (do not weaken the
   co-resident check — that is the whole point of the test).

## Acceptance

```
uv run python -m unittest tests.unit.test_zone_a_p_genenforce -v 2>&1 | tail -6
uv run python -m unittest discover tests/unit 2>&1 | tail -3
uv run python -m unittest discover tests 2>&1 | tail -3
```
Expect: `test_zone_a_p_genenforce` fully OK; unit suite OK with
**skipped=3, expected failures=1** (`test_haptic_face_golden_match` remains the
only xfail); root suite OK. Also re-confirm TASK-31's wins still hold:
```
uv run python -m unittest tests.unit.test_stereo_roundtrip_diagnostics -v 2>&1 | grep -iE "p_stereocenter_roundtrip|OK|FAIL"
```
`test_p_stereocenter_roundtrip` passes (not xfail).

## Constraints / DO NOT

- Keep the TASK-31 guard fix. Touch primarily
  `tests/unit/test_zone_a_p_genenforce.py`; do not modify `src/` enforcement
  logic to make a test pass.
- Do NOT chase the `[P@]`-on-SPL "could not be enforced" asymmetry here — just
  avoid it in fixtures and leave the note in NOTES.md for a future task.
- Run `uv run ruff format` on any file you touch. Do NOT git commit — leave
  everything (TASK-31 fix + these test changes) unstaged for review.

## On completion

Set `Status: DONE`, append a dated Log entry to `spec/worklog/NOTES.md`
(fixtures chosen, per-test retarget, acceptance results), flip TASK-31 Status to
DONE (its blocker is cleared), and update the Session-state block: suite now
skipped=3 / expected failures=1, both TASK-31 and TASK-32 done.
