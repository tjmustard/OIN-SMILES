# TASK-20: Phase 2 diagnostic — genuine P/N-stereocenter round-trip

Status: DONE — real gap confirmed (encoding-side, collapses into Phase 4). See spec/worklog/NOTES.md Log entry for 2026-07-03.
Depends on: TASK-10 (Phase 0), Stereo Phase 1 (winding plumbing, committed 6820d3a)
Suggested model: Sonnet (diagnostic); fixture step is HUMAN (chemistry)

## Goal

Determine empirically whether the OIN→XYZ generation path preserves a
stereocenter located **on the P or N atom itself** (not on a backbone carbon).
Phase 0 proved carbon `@/@@` survives; it also proved the existing
BDPP/BDNN fixtures do NOT test P/N centers. This task closes that gap and
decides whether Phase 2 needs a code fix at all.

## Prerequisite (HUMAN — chemistry, do this OUTSIDE a coding session)

Provide a new fixture at `tests/fixtures/<name>.xyz` that is a valid 3D TMC in
which a **P or N atom is a genuine CIP stereocenter** (three *distinct*
substituents + the metal bond). It must be built INDEPENDENTLY of the OIN
pipeline (otherwise the round-trip is circular). Acceptable sources:
- A published crystal structure (preferred — trustworthy geometry), or
- RDKit standalone embed: write the ligand/complex SMILES with explicit `@`/`@@`
  on the P (or N), `AddHs` + `EmbedMolecule` (ETKDGv3) + `MMFFOptimizeMolecule`,
  export XYZ. This is independent of `oinsmiles` code, so it is a fair oracle.

Recommended candidates (P is cleaner than N — P doesn't pyramidal-invert):
- **PAMP-type, monodentate (simplest):** methyl(phenyl)(2-methoxyphenyl)
  phosphine, `C[P@](c1ccccc1)c1ccccc1OC`, on a square-planar Pd/Pt with
  Cl co-ligands — geometry class you already handle well (cf. cisplatin).
  P bears Me / Ph / o-anisyl / metal → genuine stereocenter.
- **DIPAMP-type, bidentate:** (R,R)-1,2-bis[(2-methoxyphenyl)(phenyl)-
  phosphino]ethane — the classic P-stereogenic chelate; heavier, tests the
  bidentate placement path too.
Record in the fixture header comment: the molecule name, the intended CIP
label at the P/N center, and how the XYZ was produced.

Also add the human-reviewed expected OIN string as a golden:
`tests/candidate_outputs/<name>_oin.txt` (generate once with
`XYZToSMILES().convert(path)`, eyeball that the P/N carries `@`/`@@`, commit it).

## Diagnostic tests (Sonnet, once the fixture exists)

Add to `tests/unit/test_stereo_roundtrip_diagnostics.py` (same file/pattern as
TASK-10), two tests, initially `@unittest.expectedFailure`:

1. `test_p_stereocenter_roundtrip` — XYZ → OIN(1) → `OIN3DGenerator.generate()`
   → temp XYZ → OIN(2). First assert OIN(1) shows `@`/`@@` ON THE P ATOM
   (sanity; if not, the fixture or XYZ→OIN encoding is the problem, report
   that). Then assert OIN(1) == OIN(2).
2. `test_p_stereocenter_flip_inverts_cip` — take the fixture's OIN, make a twin
   with the P `@`↔`@@` flipped, generate 3D from both, and assert the two
   products have OPPOSITE CIP at the P atom via
   `Chem.AssignStereochemistryFromMol`/`AssignStereochemistryFrom3D` on
   `GeneratedStructure.mol` (the RDKit oracle, per roadmap H-1). If `.mol` is
   None for this complex, fall back to re-encoding both to OIN and asserting
   the P tags differ.

Follow TASK-10's protocol: if a test UNEXPECTEDLY PASSES, convert it to a plain
hard-assert test and record it — that means generation already preserves P/N
stereo and Phase 2 needs only the lock-in test, no code fix.

## Acceptance

```
uv run python -m unittest tests.unit.test_stereo_roundtrip_diagnostics -v
uv run python -m unittest discover tests/unit 2>&1 | tail -3   # suite still OK
```

## Outcome → next step (this is the decision point)

- **Both pass** → generation preserves P/N stereo. Phase 2 = done (lock-in
  tests only). Skip straight to Phase 3.
- **Either fails** → real gap. NOW author the fix via the HACF chain
  (`/hyper-architect` … using this diagnostic as the acceptance test). Likely
  fix sites per roadmap: fragment-SMILES→mol boundary in `_stitch_fragment`
  (`molassembler_adapter.py`), or post-embed check-and-reflect of the P/N
  stereocenter.

**Actual outcome (2026-07-03):** Both tests fail, but NOT at the round-trip
assertion — they fail at the SANITY assertion, before generation ever runs.
`XYZToSMILES().convert()` on `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` (both P
atoms genuine CIP stereocentres) produces no `@`/`@@` on either P atom at
all. Root cause: `ChiralityRecoveryUtility.recover()`
(`src/oinsmiles/core/chirality.py:155-158`) unconditionally clears the
chiral tag on any P/N atom with `total_degree < 4` in the
post-fragmentation ligand mol ("Zone A"). Since OIN's ligand-fragment SMILES
excludes the metal by construction, a directly metal-bound P/N atom always
has exactly 3 fragment-local neighbours — it is *always* Zone A, for *any*
fixture. This is not a fixture problem; it is an **XYZ→OIN encoding gap**,
not the OIN→XYZ generation gap the roadmap's Phase 2 was scoped to test.
`PseudoAtomStrategy` (`chirality.py:22`), designed to backfill the missing
4th substituent for exactly this case, is never invoked anywhere in `src/`
(dead code).

**Decision:** Phase 2 as originally scoped cannot be validated — there is
nothing for the OIN→XYZ generation path to preserve, because the signal is
already gone at encoding time. Phase 2 folds into Phase 4 (Zone-A P/N
stereo + builder decision): Phase 4 must first decide how (or whether)
Zone-A stereocenters get encoded before a Phase-2-style ETKDG experiment has
anything to test. `ROADMAP-stereo.md` Phase 2 section updated accordingly.
No code fix authored in this task (out of scope — measurement only).

## Constraints / DO NOT

- Do NOT modify `src/` in this task — it only measures.
- Do NOT build the fixture from the OIN pipeline's own output (circular).

## On completion

Set `Status:` (DONE or the outcome), append a dated Log entry to
`spec/worklog/NOTES.md` with the per-test result and the Phase-2 decision it
implies, update `ROADMAP-stereo.md` Phase 2.
