# Known limitations

Failure modes that are understood, reproducible, and *not* bugs in the component
where they surface. Each entry names the layer that would have to change.

---

## Porphyrinoid macrocycles (`macrocycle_perception`) -- mostly resolved

**History.** A metalloporphyrin's input encode carries an aromatic (or partly
localized) macrocycle, while the re-encode of its generated 3D structure localized
the whole macrocycle with E/Z markers:

```
smiles_1  [Ni_SPL].Cc1c2n{0}c(c(C)c3ccc(n{1}3)c(C)c3n{2}c(c(C)c4ccc1n{3}4)C=C3)C=C2
smiles_2  [Ni_SPL].C/C1=C2\C=CC(=N{0}2)/C(C)=C2/C=CC(=N{1}2)/C(C)=C2/C=CC(=N{2}2)/C(C)=C2/C=CC1=N{3}2
```

Both fragments are the same molecule ignoring stereo. The earlier hypothesis was a
donor-charge failure (a free-base porphyrin is a dianion; the OIN string carries no
charge, so MetalloGen builds all four N neutral and the contract mol won't sanitize).
That framing is now moot for the round trip:

1. **Aromatic-vs-localized no longer fails.** `oin/compare.py::canonical_roundtrip_key`
   (v0.3.6.1) collapses Kekule/aromatic notation, so an aromatic `smiles_1` and a
   localized `smiles_2` that describe the same skeleton match. On current `main`
   **25/46** of the registry's `macrocycle_perception` rows already round-tripped.
2. **The E/Z-only residual was a fast-path artifact, now fixed.** The remaining
   failures differed from the input by *only* E/Z slashes on the macrocycle's meso
   `C=C`/`C=N` bridges. Those bonds are ring-locked through the metal and carry no
   free E/Z; the forward encode strips them
   (`translator.XYZToSMILES.convert` -> `_clear_chelate_locked_bond_stereo`), but the
   generator's fast re-encode (`build_contract_mol` -> `get_oin_string`) skipped that
   clear, so `AssignStereochemistryFrom3D` left spurious markers. `build_contract_mol`
   now applies the same clear. Result: **34/46 (73%)** round-trip, +9, zero
   regressions. Guard: `tests/unit/test_contract_mol_chelate_ez.py`.

**Remaining residual (12/46), owned elsewhere:**
- **Pendant-bond E/Z (2):** `VAVRAN` (an exocyclic hydrazone/azo), `XIZXAG` (an
  amidine whose generated conformer has the opposite `C=N` sense). The differing bond
  is *not* metal-ring-locked, so the clear above correctly leaves it -- this is
  bond-stereo consistency (S6) / comparator (R3), not the contract mol.
- **Reduced-porphyrin skeleton/tautomer (8):** chlorins/bacteriochlorins
  (`ETPNFL`, `KECJUA`, `PRPHZN`, ...) whose generated 3D perceives a different
  sp2/sp3 saturation pattern -- a geometry/perception issue (R2), not notation.
- **Donor-H atom count (2):** `BEJSUH`, `TIDJIZ` gain/lose H at a donor (S1 domain).

---

## Ylide and radical ligands -- bond orders round-trip; residuals are downstream

Molecules whose ligands contain phosphorus ylides (`=P`) or radical carbons
(`[CH]`) -- the former `kekulize_encode_crash` bucket -- now encode without
crashing (`AGUFEN`'s PPN+ counter-cation is perceived at its correct `q=+1`) **and
their `=P`/`[CH]` bond orders now round-trip**: the contract-mol template transfer
reproduces them and the downstream radical-fill in `OINSanitizer.generate_robust_smiles`
restores any valence-deficit radical. `NAXDOI`'s `[PH]{1}(=C1[CH][CH]...)` re-encodes
byte-identical to its input.

The former `kekulize_encode_crash` starters no longer share one failure mode; each
residual is a *different*, separately-owned defect exposed once the bond orders match:

- `NAXDOI` -- its remaining diff was E/Z on a `C=C` that bridges the eta-arene and the
  ylide P (both bind Mn), i.e. metal-ring-locked. The `build_contract_mol` chelate-locked
  E/Z clear (see the macrocycle entry) removes it; the row then surfaces a downstream
  atom-count mismatch (donor-H, S1), no longer a bond-stereo failure.
- `FAMFUV` -- an sp3 atom-stereo diff (`[C@@H]` on the `[CH]`-radical ring). The forward
  encode deliberately does not enforce sp3 handedness the generator cannot reproduce
  (`translator.XYZToSMILES.convert`); the fast re-encode's `AssignStereochemistryFrom3D`
  does. This is the atom-stereo analog of the E/Z fix and belongs to bond/atom-stereo (S6),
  not the contract-mol E/Z path.
- `DAVQIA` -- conformer generation / atom count (R2 / S1).

Encoder-side and bond-order-side this is closed. The remaining work is atom-stereo (S6),
donor-H atom count (S1), and conformer robustness (R2).

---

## Borane and carborane clusters

`get_lig_mol` cannot perceive bond orders for electron-deficient cages (`OSENOR`,
and the `carborane_unsupported` bucket). It already fails with a specific message
(`get_lig_mol failed for ligand fragment #4 (SMILES: 'B')`), and the charge sweep
does not rescue it: the multi-centre bonding these cages need is outside the
two-centre model `AC2mol` implements.
