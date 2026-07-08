# S4 — Eta-winding canonicalization for automorphic rings + eta-slot placement

Branch: `feature/roundtrip-eta-winding` · Read `docs/handoffs/README.md` first.

## Mission

Round trips fail on winding markers (`{n>}` vs `{n<}`) for eta rings whose
substituents are all identical (Cp*, benzene, C6Me6) — 80 registry rows flag a
winding difference. For a fully automorphic ring, CW vs CCW traversal describes
the SAME structure (the ring has in-plane C₂ axes), so this is notation
ambiguity, not a wrong diastereomer — but winding IS load-bearing for
substituted rings (TiCat3/4 rac/meso, the v0.3.4 eta-winding fix). Success =
automorphic-ring complexes round-trip regardless of embed orientation, while
every substituted-ring winding guard stays green.

## Evidence pack

`CAHKEE_comp_0` [Rh_TET] — Cp*Rh(carbene)Cl₂, ONLY difference is one marker:

```
1: ...Cc{1}1c{1>}(C)c{1}(C)c{1}(C)c{1}1C...
2: ...Cc{1}1c{1<}(C)c{1}(C)c{1}(C)c{1}1C...
```

`SOJMIQ_comp_0` [Co_TPL] — second gap, same family: a BPh₄⁻-like borate with
one η-bound phenyl; input and re-encode put the eta slot markers on DIFFERENT
(chemically equivalent) phenyls of the fragment, so the strings differ although
the structure is the same:

```
1: c1ccc(B(c{2}2[cH]{2}[cH]{2>}[cH]{2}[cH]{2}[cH]{2}2)(c2ccccc2)c2ccccc2)cc1
2: c1ccc(B(c2ccccc2)(c2ccccc2)c{2}2[cH]{2}[cH]{2}[cH]{2>}[cH]{2}[cH]{2}2)cc1
```

Fresh rows also flagged winding+skeleton (`RAJNOI`, `TIKRUC`, `MEFWAY`) — the
skeleton part is S2's diene bug; re-check them after S2 lands. Many of the 80
are stale-vintage; use the post-Phase-0 `CASE_REGISTRY.md` list.

## Where the logic lives

- Encoder winding emission: `src/oinsmiles/utils/oin_aligner.py` (per-eta-slot
  winding via the ring's actual metal→centroid axis — the 2026-07-06 fix; RC1/
  RC2 eta-ring canonicalization: scoped rank swap by canonical ring SMILES +
  heading atom = lowest `Chem.CanonicalRankAtoms`).
- Comparator: `src/oinsmiles/oin/compare.py` — `winding_canonical_key` already
  collapses which-of-two-equivalent-RINGS-is-slot-0 relabeling via a winding
  MULTISET, but a single automorphic ring flipping `>`↔`<` changes the multiset
  and fails. Its docstring explicitly documents the rac/meso and enantiomer
  cases that MUST keep failing.

## Design guidance (verify, then choose the layer)

Preferred: fix at the ENCODER (canonical tie-break), not the comparator — the
OIN string should be deterministic for the same structure. For a ring whose
eta atoms are all in one automorphism class (check with
`Chem.CanonicalRankAtoms(breakTies=False)` on the ligand fragment): the winding
direction is not structure-bearing, so always emit a fixed direction (e.g. `>`).
For substituted rings, keep the geometric winding exactly as today.
A comparator-side collapse gated on the same automorphism test is acceptable as
a belt-and-braces second layer, but must NOT collapse substituted rings.

SOJMIQ: extend the canonical-representative idea (see
`OINSanitizer.canonical_donor_representative`, the Track-A1 symmetric-donor
fix in `oin_aligner.py`) from single donor atoms to eta ring SETS: among
automorphic candidate rings, place markers on the canonical one.

## Verify-first steps

1. Confirm the equivalence claim on CAHKEE: generate twice with different seeds
   (or flip the marker by hand) → same 3D structure/RMSD → notation-only.
2. Find the tie-break point in `oin_aligner.py` where winding is computed;
   check what breaks ties today for automorphic rings (probably atom order =
   nondeterministic across embeds).
3. Run the existing eta-winding guards BEFORE changing anything:
   `tests/unit/` TiCat/eta tests + `test_roundtrip_canonical_key.py` — know
   your regression floor.

## Files

- **Own:** `src/oinsmiles/utils/oin_aligner.py`, `src/oinsmiles/oin/compare.py`
  (+ new test file `tests/unit/test_automorphic_ring_winding.py`).
- **Read-only:** `utils/xyz2mol.py` (S3), `metallogen_adapter.py` (S1/S2).
- **Regression floor:** every existing winding/eta test, especially the
  TiCat3/4 rac-meso guards and `tests/unit/test_roundtrip_canonical_key.py`.

## Acceptance

- `CAHKEE_comp_0` round-trips across ≥3 embed seeds.
- `SOJMIQ_comp_0` round-trips (marker placement canonical).
- New guards: automorphic ring (Cp*) winding stability + a substituted-ring
  counter-case that MUST still distinguish `<`/`>` (rac vs meso).
- Full unit suite green on BOTH rdkit 2025.09.3 and 2026.3.3 (canonical ranks
  can differ across versions — the v0.3.5 E/Z work hit exactly this).
