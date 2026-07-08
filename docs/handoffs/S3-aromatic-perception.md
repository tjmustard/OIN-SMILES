# S3 — Aromatic/quinoid garbling in the XYZ→OIN encoder (both directions)

Branch: `feature/roundtrip-aromatic-perception` · Read `docs/handoffs/README.md` first.

## Mission

The encoder's aromatic perception breaks on quinoid/macrocyclic conjugation, in
three escalating forms: (1) generated structures re-encode with garbled mixed
notation (`c1c=c(C)...` — lowercase aromatic atoms carrying explicit double
bonds), (2) porphyrin-type macrocycles fully dearomatize into localized E/Z
strings with `[CH][CH]` radical-looking carbons — sometimes already in the
INPUT encode, (3) hard crashes: `Can't kekulize mol` (13 rows) and
`xyz2mol failed: 'NoneType' object has no attribute 'GetAtoms'` (7 rows).
47 + 13 + 7 + 4 registry rows. Success = the evidence molecules round-trip (or
crash cases produce honest, specific errors), classes `garbled_aromatic`,
`kekulize_encode_crash`, `xyz2mol_none_crash` ≈ 0.

## Evidence pack

`KIYWUM_comp_0` [Au_LIN] — mesityl ring garbled ONLY in the re-encode:

```
1: Cc1cc(C)c(-c2c{0}n(...)nn2C)c(C)c1
2: Cc1c=c(C)c(-c2c{0}n(...)nn2C)=c(C)c=1     # aromatic c with explicit = bonds
```

`TIPDIG_comp_0` [Ni_SPL] — porphyrinoid: input aromatic `c1c2n{0}c(...)`;
generated re-encode fully dearomatized `/C1=C2\[CH][CH]C(=N{0}2)/...` with
bare `[CH]` carbons (radical notation) and spurious E/Z decorations.

`GOBSUL_comp_0` [Zn_SPL] — macrocycle where even smiles_1 contains
`C1[CH]C(=c2...)` — the INPUT encode is already semi-garbled, so this is not
purely a generated-structure problem.

Fresh same-class: `SOJMIQ`-adjacent rows `DEDWEP`, `HUXJOB`, `XOJHIN`,
`QIGPON`, `SEJPEE` (Rh_SPL), `HUVCUW`, `YAXVOJ` (Ti_TET), `YIQSIB` (Mn_OCT),
`XIRDAF` (Rh_OCT), `ZOTVIP` (Ni_OCT), `HOSXUJ` (Zn_OCT), `NOBYOU` (Fe_OCT).
Crashes: `NAXDOI_comp_0`, `FAMFUV_comp_0` (kekulize, fresh), `HIMQIF_comp_0`
is carborane — NOT yours (wontfix-docs). Full list: `CASE_REGISTRY.md`.

## Prior art / known context

- **Track C (v0.3.5)** fixed the quinoid PARSE crash on the generation side:
  `_dearomatize_stuck_rings` in `src/oinsmiles/generator3d/process.py` clears
  aromaticity on only the stuck ring(s) (aromatic atom with exocyclic double
  bond) so normal benzene rings keep alternating orders. Guard:
  `tests/unit/test_quinoid_ligand_parse.py`.
- Its documented residual is exactly your class: “the 2-iminopyridine ligand
  can't round-trip cleanly: contract-mol re-encode gets benzene wrong (quinoid)
  but imine `=` right; pure-geometry re-encode gets benzene right but drops the
  imine C=N” (stale cohort `ABERIK/ABEROQ/ABERUW/ACUYUU/AFAMEB/AFAMIF/AFAMOL/
  AFIXUJ` — Ti/Hf amidinate & 2-iminopyridine).
- The re-encode path: harness calls `get_oin_string(gen.mol, coords)`
  (`src/oinsmiles/utils/xyz2mol.py`) when the generator returns a bonded mol,
  else full `XYZToSMILES.convert`. The garbling appears in the first path; the
  hard crashes in the second (input side, `get_lig_mol` / kekulize during
  fragment SMILES building).

## Verify-first steps

1. Repro `KIYWUM_comp_0` (`--only`). Dump the fragment mol right before SMILES
   emission in `get_oin_string` — find where aromatic flags and bond orders
   disagree (garbled `c=c` means aromatic-flagged atoms with KEKULE double
   bonds and no matching ring perception).
2. Repro `TIPDIG_comp_0` — check whether `build_contract_mol` transferred
   template bond orders (S2 owns that function; if the root cause lands there,
   coordinate via PR notes rather than editing it — your fix may be to make the
   ENCODER robust to both inputs: re-aromatize/sanitize the fragment before
   emission).
3. For the kekulize crashes (`NAXDOI_comp_0`): the input-side ligand builder
   raises during Kekulize on 5-ring aromatic systems (`Unkekulized atoms`
   lists). Consider the same stuck-ring-scoped strategy as
   `_dearomatize_stuck_rings`, applied in the encoder's `get_lig_mol` path.
4. `xyz2mol_none_crash` (7 rows): `get_tmc_mol` returned None — triage what
   these molecules share before designing anything.

## Files

- **Own:** `src/oinsmiles/utils/xyz2mol.py`, `src/oinsmiles/generator3d/process.py`
  (+ new test file `tests/unit/test_aromatic_reencode.py`).
- **Read-only:** `metallogen_adapter.py` (S1/S2), `utils/oin_aligner.py` and
  `oin/compare.py` (S4).
- **Regression floor:** `tests/unit/test_quinoid_ligand_parse.py`.

## Acceptance

- `KIYWUM_comp_0` re-encodes with a clean aromatic mesityl ring; keys match.
- `NAXDOI_comp_0`/`FAMFUV_comp_0` either encode successfully or fail with a
  specific, documented limitation error (not a bare kekulize traceback).
- Porphyrinoids (`TIPDIG`): forward encode is stable (encode the SAME input
  twice → identical string) and the re-encode either matches or the residual
  is documented as a template-transfer issue for S2's pipeline.
- Guard tests for: quinoid ligand round-trip, garbled-aromatic re-encode,
  kekulize-crash input. Full unit suite green; stale 2-iminopyridine cohort
  spot-checked (`ABERIK_comp_0`).
