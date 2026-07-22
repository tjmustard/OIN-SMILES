# Conformer ensembles — B2 convergence fixtures

Multi-frame `crest_conformers.xyz` files: genuinely distinct 3D geometries (conformers) of
one isomer. `tests/integration/test_conformer_convergence.py` asserts every frame of an
ensemble encodes to the **same** canonical key
`winding_canonical_key(normalize_oin_for_comparison(XYZToSMILES().convert(frame)))`.

This is the *true-conformer* companion to the rigid rotation/translation guard in
`test_conformer_invariance.py` (`tests/fixtures/conformer_set/`). There the gate is
byte-identical raw strings (same geometry, just reoriented); here genuine distinct geometries
carry benign rotational slot drift, so the gate is the canonical **key**, not the raw string.

## Provenance

Frames are frame-capped copies (cap = 6, ascending energy = CREST file order) of the ensembles
produced by the CREST sweep in the sibling worktree
`oin-conformer-invariance/conformer_crest_sweep_gxtb_reopt/crest_<MOL>/crest_conformers.xyz`
(`crest --noreftopo --mquick`, `gfnff` search, `gxtb` preopt/reopt). That sweep directory is a
gitignored artifact; these curated, size-bounded copies are committed here. Verdicts are from
its `conformer_invariance_report.json`. Regenerate with
`scratchpad/build_ensemble_fixtures.py` (see `manifest.json:generated_with`).

## Selection

- **Fast tier** (always run): `CisPlatin` (Pt, 4 frames), `OPAGES_comp_0` (Hg, 5 frames) —
  small, exercise a couple of genuinely distinct geometries quickly.
- **Full tier** (`OIN_CONVERGENCE_FULL=1`): `QEXRAP` (Pd), `XIXPEB` (Os), `MAZJED` (Ti),
  `COKGAN` (Ir), `NODNOL` (Cr), `KARQUQ` (V), `QUSJAQ` (Fe), all `invariant`; plus:
  - **`BEPCAC_comp_0`** (Ni) — a *regression guard for the B1 electronic geometry prior*.
    Without B1 this d8 complex splits SPL/TET across conformers; with B1 it converges. Revert
    B1 and this fixture fails.
  - **`CETDAI_comp_0`** (Hf), verdict `notation-drift` — its committed frames span **2 distinct
    raw OIN encodings** (symmetric guanidinate donor slot order) that share **one** canonical
    key. The convergence test asserts exactly that (`>=2` raw, `1` key): a positive demonstration
    that the canonical key absorbs benign notation drift.

## Excluded

- `YIJJEH_comp_0` (Ni, review-divergent): 3 distinct OIN → **2 distinct keys** (a genuine
  square-planar/tetrahedral geometry ambiguity, τ4 ≈ 0.5). Not a convergence fixture — it is a
  real divergence that would (correctly) fail the convergence gate.
- `single-conformer` molecules (e.g. `TransPlatin`): CREST produced < 2 conformers, so there is
  nothing to converge.

## Format

Standard concatenated multi-frame XYZ: `N` / comment (CREST energy, Hartree) / `N` atom lines,
repeated per frame, no blank separators. Frames are ordered by ascending energy.
