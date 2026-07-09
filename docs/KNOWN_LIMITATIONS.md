# Known limitations

Failure modes that are understood, reproducible, and *not* bugs in the component
where they surface. Each entry names the layer that would have to change.

---

## Porphyrinoid macrocycles re-encode localized (`macrocycle_perception`)

**Symptom.** A metalloporphyrin's input encode carries an aromatic macrocycle,
while the re-encode of its generated 3D structure is fully localized with E/Z
markers, so the round-trip key does not match:

```
smiles_1  [Ni_SPL].Cc1c2n{0}c(c(C)c3ccc(n{1}3)c(C)c3n{2}c(c(C)c4ccc1n{3}4)C=C3)C=C2
smiles_2  [Ni_SPL].C/C1=C2\C=CC(=N{0}2)/C(C)=C2/C=CC(=N{1}2)/C(C)=C2/C=CC(=N{2}2)/C(C)=C2/C=CC1=N{3}2
```

Both fragments are the same molecule ignoring stereo. RDKit re-parses *neither* as
aromatic; the difference is only which bonds got E/Z labels.

**This is not an aromatic-perception bug in the encoder.** A free-base porphyrin is
a dianion. `get_tmc_mol` perceives that: in `tmc_mol` two of the four pyrrole
nitrogens carry a `-1` formal charge, and because the encoder sanitizes with the
metal still attached through its dative bonds, RDKit sees the 18-pi macrocycle and
marks it aromatic. The ligand fragments inherit those flags.

An OIN string carries no formal charges, so `n{0}` is ambiguous between a neutral
pyridine-type nitrogen and an anionic pyrrolide one (the same ambiguity that makes
a bare `N{n}` mean both nitride and ammine). MetalloGen therefore builds all four
nitrogens neutral. The resulting contract mol cannot even sanitize:

```
Explicit valence for atom # 7 N, 4, is greater than permitted
```

With no valid valence model there is no aromaticity to perceive, so the re-encode
localizes the macrocycle. Verified on `BEGLUU_comp_0`:

```
tmc_mol       Ni nbrs = [(N,-1,DATIVE), (N,0,DATIVE), (N,-1,DATIVE), (N,0,DATIVE)]
contract mol  Ni nbrs = [(N, 0,DATIVE), (N,0,DATIVE), (N, 0,DATIVE), (N,0,DATIVE)]
```

**Forward-encode stability holds.** Encoding the same input three times yields the
identical string, so nothing here is non-deterministic.

**What would fix it.** The donor's charge has to survive the OIN round trip -- either
the format marks an anionic donor, or the OIN -> m-SMILES conversion infers it from
the metal's oxidation state. That is the donor-charge layer, not the encoder. Both
`_perceive_aromaticity_if_absent`-style re-perception in `get_oin_string` and a
whole-mol `Chem.SetAromaticity` were tried and cannot work: the contract mol fails
to sanitize before any aromatic model runs.

---

## Ylide and radical ligands encode but do not round-trip

Molecules whose ligands contain phosphorus ylides (`=P`) or radical carbons
(`[CH]`) -- the former `kekulize_encode_crash` bucket -- now encode without
crashing, and `AGUFEN`'s PPN+ counter-cation is perceived at its correct `q=+1`.
They still fail the round trip downstream: MetalloGen cannot reproduce these
bond orders, so the comparison fails on atom count or string.

Encoder-side this is closed. The remaining work is in 3D generation.

---

## Borane and carborane clusters

`get_lig_mol` cannot perceive bond orders for electron-deficient cages (`OSENOR`,
and the `carborane_unsupported` bucket). It already fails with a specific message
(`get_lig_mol failed for ligand fragment #4 (SMILES: 'B')`), and the charge sweep
does not rescue it: the multi-centre bonding these cages need is outside the
two-centre model `AC2mol` implements.
