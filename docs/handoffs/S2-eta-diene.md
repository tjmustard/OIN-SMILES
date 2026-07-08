# S2 — η²-alkene/diene (COD) double-bond localization in generated ligands

Branch: `feature/roundtrip-eta-diene` · Read `docs/handoffs/README.md` first.

## Mission

Complexes with η²-bound alkenes — above all **1,5-cyclooctadiene (COD) on
Rh/Ir**, ubiquitous in tmCAT — re-encode with the C=C bonds moved OFF the
metal-bound carbons onto the ring backbone. 40 registry rows; Rh's overall 66%
failure rate is largely this class. Success = COD complexes round-trip
byte-stable on the eta fragment; `eta_diene_localization` ≈ 0 in the post-fix
registry.

## Evidence pack

`GASBIN_comp_0` [Rh_SPL] — input vs generated re-encode:

```
1: ...[CH]{1}1=[CH]{1>}CC[CH]{2>}=[CH]{2}CC1...      # η2,η2-COD: C=C on the bound carbons
2: ...[CH]{1}1[CH]{1>}[CH2]=[CH2][CH]{2>}[CH]{2}[CH2]=[CH2]1...   # C=C migrated to the CH2-CH2 backbone!
```

`PENGAT_comp_0` [Ir_TBP] — same signature. Also fresh: `MEFWAY_comp_0`
[Ru_TET], `XULBAJ_comp_0` [Ru_OCT], `ZOXLIJ_comp_0` [Ru_TBP], `RAJNOI_comp_0`
[Rh_TET], `TIKRUC_comp_0` [Os_TET] (the last three also flag a winding
difference — the winding side belongs to S4; your job is the bond orders).
Full list: `CASE_REGISTRY.md` → `eta_diene_localization`. Smoking gun:
generated string contains `[CH2]=[CH2]` inside a ring where the input has
`[CH]{n}=[CH]{n}`.

## Root cause — HYPOTHESIS (verify first; this class has burned two handoffs)

The η³-allyl "double-bond loss" (fixed in v0.3.5) had this mechanism: an OIN
ligand atom bound to the (stripped) metal is under-valent → RDKit gives the
template atom a RADICAL electron → `_flatten_template`'s flattened pattern
keeps that radical → `GetSubstructMatch` treats radical count as a match
constraint → the match silently fails → `build_contract_mol` transfers NO bond
orders → re-perception localizes bonds wrong. Fix was clearing
radicals/NoImplicit/ExplicitHs in `_flatten_template`
(`src/oinsmiles/generation/metallogen_adapter.py`), guarded by
`tests/unit/test_contract_mol_allyl_transfer.py`.

η²-diene carbons are `[CH]` with an explicit H and a double bond — plausibly a
*different* failure in the same transfer pipeline (e.g. the template match
succeeds but maps symmetric ring atoms to the wrong automorph, so bond orders
transfer onto rotated positions — COD has a 4-fold-ish symmetric skeleton; or
the match fails on H-count rather than radicals). **Do not assume; measure:**

1. Repro `GASBIN_comp_0` with `--only`; capture `gen_result.mol`.
2. Instrument `build_contract_mol`: does `GetSubstructMatch(flattened_template)`
   return a mapping? If empty → transfer skipped (allyl-style). If non-empty →
   print the atom mapping and check whether template C=C pairs land on the
   generated ring's CH2–CH2 positions (automorphism mis-map).
3. Depending on 2: either extend the `_flatten_template` clearing, or make the
   match/transfer stereo-aware of WHICH automorph to pick (e.g. anchor the map
   on the binding-slot atoms — they are known from the OIN vectors).

## Files

- **Own:** `src/oinsmiles/generation/metallogen_adapter.py` — ONLY
  `_flatten_template`, `build_contract_mol`, `_oin_fragment_templates` (+ new
  test file `tests/unit/test_contract_mol_diene_transfer.py`).
- **Read-only:** `convert_parsed_to_msmiles` (S1 owns it), `utils/xyz2mol.py`
  (S3), `utils/oin_aligner.py` / `oin/compare.py` (S4).
- **Regression floor:** `tests/unit/test_contract_mol_allyl_transfer.py` (4
  tests) must stay green — the allyl fix lives in the code you are editing.

## Acceptance

- `GASBIN_comp_0`, `PENGAT_comp_0` round-trip with canonical keys matching and
  the eta fragment's `=` on the bound carbons.
- New guard test: COD-like η²,η² template transfer (mirror the allyl test
  structure).
- Regression: allyl cases `ABAZEK_comp_0`, `ACALOI_comp_0` still round-trip
  (they passed after v0.3.5); full unit suite green.
- Cases that fail for a DIFFERENT residual reason after your fix (e.g. RMSD-996
  mapping on `ABETIK`-like bulky systems) — document in the PR body and leave
  to S5.
