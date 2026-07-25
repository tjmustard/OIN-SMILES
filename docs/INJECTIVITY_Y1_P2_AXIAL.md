# Injectivity Y1 · P2 — Axial / atropisomeric chirality (ENCODER-BLIND)

**Status:** measurement only. Verdict: **CONFIRMED encoder-blind (total)**. Parent:
`INJECTIVITY_Y1_OVERVIEW.md`.

## The isomer pair

`PdCl₂(BINAP)` — the ligand is 2,2′-bis(diphenylphosphino)-1,1′-binaphthyl, whose stereochemistry
is a hindered biaryl **atropisomeric axis** (`R`/`S`, aka M/P). The z-mirror of the R-BINAP
fixture is the S-BINAP complex. The independent geometric oracle confirms distinctness: min
proper-rotation mirror RMSD **4.02 Å**.

## Measurement

`convert(R-BINAP)` and `convert(mirror)` are **byte-identical**:

```
[Pd_SPL].c1ccc(P{0}(c2ccccc2)c2ccc3ccccc3c2-c2c(P{1}(c3ccccc3)c3ccccc3)ccc3ccccc23)cc1.[Cl]{2}.[Cl]{3}
```

`raw_equal = True`, `key_equal = True` → **ENCODER-BLIND (total)**: not even the raw string
separates the two atropisomers.

## Mechanism

The forward encode perceives only C=C/C=N double-bond E/Z (`Chem.DetectBondStereochemistry`,
`core/translator.py:115`) — there is **no single-bond atropisomeric-axis perception anywhere**
in `core/translator.py` or `core/chirality.py`, so the biaryl `-` bond carries no `/`/`\`
descriptor. Notably the independent RDKit `FindPotentialStereo` fingerprint on the metal-free
ligand is also empty here: standard cheminformatics does not perceive this axis either, which is
part of why the encoder does not.

## Verdict & disposition

**CLOSED as an opt-in, round-tripping fix (Wave 2).** NB the mechanism differs from the plan's
assumption: RDKit does *not* perceive this axis from pure 3D (per the RDKit Book an atropisomer
bond is only marked when a neighbour bond is *wedged*; `FindPotentialStereo` returns nothing on
the BINAP ligand). The recoverable signal is the **signed biaryl dihedral**, computed in
`src/oinsmiles/oin/axial.py`.

Behind `OIN_EMIT_AXIAL` (default OFF → output byte-identical):

- the encoder appends a **canonical** axial token — invariant under atom renumbering and proper
  rotation, flipping only under reflection — so R-BINAP emits `|ax:-|` and S-BINAP `|ax:+|`;
- a **stereogenicity gate** suppresses axes whose ring end has symmetry-equivalent ortho
  neighbours (locally C2 ⇒ achiral however twisted). Without it the encoder would be
  *over-sensitive*, claiming chirality that does not exist; the independent geometric oracle
  agrees with the gate exactly (mirror RMSD 0.000 Å on such cases);
- the round-trip **key folds** the token, so the batch harness is unaffected either way;
- the **generator honours it**: an axial-aware pass in `_select_by_geometry` plus an axial-aware
  acceptance predicate. Measured A/B (`tools/injectivity/axial_roundtrip_ab.py`): **2/2 with the
  pass vs 1/2 without** — the baseline returns the same handedness whatever is requested.

So the axis now survives `XYZ → OIN → 3D` and emitting it no longer trades a false positive for a
false negative. It remains **opt-in** pending a broader A/B than the single fixture pair (n=2).
See `docs/INJECTIVITY_Y2_FEASIBILITY.md`. Guards: `tests/unit/test_axial_emit.py` (emit +
canonicality + over-sensitivity), `tests/integration/test_axial_roundtrip.py` (round trip),
`tests/unit/test_injectivity_probes.py::test_axial_is_encoder_blind` (default-off blindness).
Reproduce: `PYTHONPATH=$PWD/src python -m tools.injectivity.twin_collision tests/fixtures/PdCl2-R-BINAP.xyz`.
