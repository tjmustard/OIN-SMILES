# Injectivity Y1 · P1 — Metal-centre Δ/Λ chirality (KEY-BLIND)

**Status:** measurement only. Verdict: **CONFIRMED key-blind** (documented-deferred limitation,
now demonstrated through `convert()` on real geometry). Parent: `INJECTIVITY_Y1_OVERVIEW.md`.

## The isomer pair

`fac-Ir(ppy)₃` is a C₃-symmetric tris(bidentate) octahedral complex — a textbook Δ/Λ chiral
species with **achiral ligands**, so its only stereochemistry is the metal-centred handedness.
Its z-mirror is the opposite (Λ vs Δ) enantiomer. The independent geometric oracle confirms
distinctness: min proper-rotation mirror RMSD **3.19 Å** (achiral controls sit at ~0.05 Å).

## Measurement

| | base (Δ) | mirror (Λ) |
|---|---|---|
| raw OIN | `[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}...n{1}1.c{2}...n{4}1` | `[Ir_OCT].c{5}...n{0}1.c{1}...n{3}1.c{2}...n{4}1` |

`raw_equal = False`, but `key_equal = True` → **KEY-BLIND**. The raw strings differ *only* by
slot renumbering (which of the three identical ppy ligands got which absolute slot), which the
round-trip key deliberately folds via `_polyhedron_signature` (`oin/compare.py`). The batch
harness gates on the key, so two genuine Δ/Λ enantiomers **round-trip as identical** — a
false positive in the dangerous cell.

## Mechanism

- The encoder emits **no** metal stereo descriptor: `OINInlineHandler.generate_inline_string`
  builds the metal token as bare `[Metal_GEO]` (`oin/inline.py:110`); the `@…` group in
  `METAL_REGEX` (`oin/inline.py:51`) is a vestigial reader with no producer.
- The key strips any `@SP/@OH/@TB` label (`_METAL_STEREO_RE`, `oin/compare.py:55,91`) and folds
  slot relabeling (`_polyhedron_signature`, `:399`) — deferred pending a *reproducible*
  encoder-side metal stereo descriptor (`compare.py:61-66,453`).

This corrects the false confidence in `tests/integration/test_isomer_divergence.py::`
`test_metal_stereo_raw_only`, which asserts metal-@ "survives in the raw string" using
**hand-fed** `@SP1`/`@SP2` synthetic OINs the encoder never actually produces.

## Verdict & disposition

**Recoverable, deferred to v0.4.5 (NOT permanent).** Wave 2 (`docs/INJECTIVITY_Y2_FEASIBILITY.md`)
showed RDKit's `AssignStereochemistryFrom3D` *does* perceive the metal Δ/Λ configuration from the
3D coordinates — the octahedral `_chiralPermutation` flips between enantiomers (fac-Ir(ppy)₃:
10 vs 8). So the information is recoverable; what is missing is a **canonical, orientation-invariant
donor ordering** to turn that raw permutation into a reproducible `@OHn`/`@SPn` token (RDKit lists
non-tetrahedral CIP/canonicalization as "totally missing"), plus a generator that can set it. That
canonical-ordering problem is exactly the v0.4.5 canonical-string work, so P1 is filed there rather
than as a permanent limitation. Guard:
`tests/unit/test_injectivity_probes.py::test_metal_delta_lambda_is_key_blind`
(+ aspirational `test_metal_chirality_should_diverge_at_key`). Reproduce:
`PYTHONPATH=$PWD/src python -m tools.injectivity.twin_collision tests/fixtures/fac-Ir(ppy)3.xyz`.
