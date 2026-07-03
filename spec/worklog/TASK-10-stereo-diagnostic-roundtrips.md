# TASK-10: Stereo diagnostic round-trip tests (expected failures)

Status: DONE
Depends on: TASK-04 (expected strings must be in v3.7 descriptor-free style)
Suggested model: Sonnet

## Goal

Add three round-trip tests that EXPOSE the currently-silent stereo loss on the
OIN→XYZ (generation) side. They are marked `@unittest.expectedFailure` so the
suite stays green while documenting the gap; they become the acceptance tests
for ROADMAP-stereo.md Phases 1–4 (each phase flips one to a hard assert).

## Context (no prior repo knowledge needed)

The XYZ→OIN direction encodes ligand P/N stereocenters as `@/@@` in the
fragment SMILES and eta-ring winding as `{n>}`/`{n<}` slot markers. The
OIN→XYZ direction (`OIN3DGenerator` in `src/oinsmiles/generation/engine.py`)
currently drops these: winding is not captured by `SLOT_REGEX`
(`src/oinsmiles/oin/inline.py:44`), P/N CIP codes are never read during
generation, and haptic ring FACE is chosen only by flipping the ring normal
toward the metal (`src/oinsmiles/generation/molassembler_adapter.py:627`).
So a stereo-bearing OIN string may generate a structure that re-encodes to a
DIFFERENT string. These tests measure exactly that.

Available fixtures (all round-trip-verified for geometry in v0.2.1):
- `tests/fixtures/PdCl2-RR-BDPP.xyz` — P-stereocenters (chiral phosphine)
- `tests/fixtures/PdCl2-RR-BDNN.xyz` — N-stereocenters (chiral amine)
- `tests/fixtures/ferrocene.xyz` — η5-Cp sandwich
Golden OIN strings (post-TASK-04): `tests/candidate_outputs/bdpp_oin.txt`,
`bdnn_oin.txt`, `ferrocene_oin.txt`.

Public API (check `src/oinsmiles/__init__.py` for exact signatures):
- XYZ→OIN: `XYZToSMILES` … `.convert()` → OIN string
- OIN→XYZ: `OIN3DGenerator` … `.generate(oin)` → `GeneratedStructure` with
  `.xyz` (XYZ block string) and `.mol` (RDKit mol or None).

## Files to touch

- NEW: `tests/unit/test_stereo_roundtrip_diagnostics.py` (only file)

## Steps

Create the new test file with a class `TestStereoRoundTripDiagnostics`
containing three tests, each `@unittest.expectedFailure` and each with a
docstring naming the roadmap phase that will fix it:

1. `test_chiral_p_roundtrip` (fixed by Phases 1/2/4):
   XYZ (`PdCl2-RR-BDPP.xyz`) → OIN (assert it contains `@` ligand tags —
   sanity, should pass today) → `OIN3DGenerator.generate()` → write the
   `.xyz` to a temp file → XYZ → OIN again → assert the second OIN string
   equals the first (including `@/@@` tags).
2. `test_chiral_n_roundtrip`: same with `PdCl2-RR-BDNN.xyz`.
3. `test_haptic_face_winding` (fixed by Phase 3): take the ferrocene golden
   OIN (`tests/candidate_outputs/ferrocene_oin.txt`), and construct a second
   string with one ring's winding flipped (`{0>}` → `{0<}`). Generate 3D from
   both, re-encode both to OIN, assert the two re-encoded strings differ in
   that ring's winding marker. (Plain ferrocene's faces are homotopic, so
   also acceptable per implementer judgment: assert simply that winding
   markers survive the round trip at all — today they do not.)

Practical notes:
- Generation can take ~10-60 s per complex (embedding + fallbacks); keep the
  test count at 3.
- If `OIN3DGenerator.generate()` raises for a fixture (not just produces
  wrong stereo), catch it and `self.fail(...)` with the error — a crash is
  also a diagnostic result; `expectedFailure` will record it.
- Use `tempfile` for intermediate XYZ files; do not write into the repo.

## Acceptance (exact commands)

```
uv run python -m unittest tests.unit.test_stereo_roundtrip_diagnostics -v
```
Expected: 3 tests, each reported as **expected failure** (or unexpectedly
passing — if one PASSES, report it in NOTES.md: it means that pathway already
preserves stereo and the roadmap phase can be downgraded to "add hard test").

```
uv run python -m unittest discover tests/unit 2>&1 | tail -3
```
Expected: suite still green (expected failures don't fail the run).

## Constraints / DO NOT

- Do NOT modify any source under `src/` — this task only measures.
- Do NOT mark tests as `skip`; use `expectedFailure` so regressions/progress
  are visible in every run.

## Out of scope

- Fixing the stereo loss (ROADMAP-stereo.md Phases 1–4, each via its own
  HACF MiniPRD).

## On completion

Set `Status: DONE`, append a Log entry to `spec/worklog/NOTES.md` recording
which of the three failed vs unexpectedly passed and with what mismatch —
that data drives the Phase 1–4 MiniPRDs.
