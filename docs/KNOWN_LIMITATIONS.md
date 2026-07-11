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

---

## Some complexes yield no conformer (`no_conformers`)

The v0.3.7 R2 triage ran all 36 rows the post-S6 registry filed under `no_conformers`
serially (the registry over-counts this class: a *concurrent* sweep fabricates the
error, so several rows are contention flakes that generate on a serial retry). After
the two R2 fixes below, **27 of the 36 generate a conformer**; the remaining **9 are
genuine and unfixable in the generator's owned layers.** They fall into three groups,
each of which surfaces where the vendored MetalloGen bond-order/valence model cannot
build a sane molecule for distance geometry.

**Group 1 -- neutral L-donor over-valence (the donor-charge gap above, at embed time).**
`FUVNER`, `GEZKAZ`, `VIBRIK`. An OIN string carries no formal charges, so a neutral
two-electron donor (amine N, aqua O, triflate O) bonded to the metal gains one bond
beyond its neutral valence and the embed mol will not sanitize:

```
GEZKAZ_comp_1  [Zn_TET].[Cl]{0}.[Cl]{1}.[Cl]{2}.[OH2]{3}
               Explicit valence for atom # 4 O, 3, is greater than permitted
VIBRIK_comp_0  [Fe_TET].CN{0}(C)CC(C)(C)N{1}(C)C.[Cl]{2}.[Cl]{3}
               Explicit valence for atom # 5 N, 4, is greater than permitted
```

This is the same root cause as the porphyrin entry -- the donor's charge/dative
character has to survive the OIN round trip -- and would be fixed in the OIN ->
m-SMILES donor-charge layer, not in the embed.

**Group 2 -- exotic bond orders the two-centre model can't perceive.** `DAHXOB`
(hypervalent-S ylide `C=C=S=C`), `MEDDUV` (cyclophosphazene `P=N` ring), `IREPAX`
(diphosphene `P=P`), `DOFCAE` (aromatic boron `[b]`), and `HURGOS` (a fused aromatic
that will not kekulize). PuLP / `get_valid_molecule` returns an invalid or `None`
perception, so `ff_clean` raises deep in the FF setup (`cannot unpack non-iterable
bool object`, `'NoneType' object is not subscriptable`, `Can't kekulize mol`). These
are the same class as the ylide/radical and carborane limitations -- a perception
gap, not an embed-parameter one.

**Group 3 -- geometry the builder cannot realize.** `BOBJIM` (`Sc_OCT`): the embed
succeeds but every FF-cleaned geometry fails the ligand-collision / adjacency
validity checks, so no conformer survives cleanup. A genuine geometry-realization
gap for that ligand set.

### What R2 *did* fix

* **`_apply_double_bond_stereo` no longer forces an over-valent double bond**
  (`generator3d/embed.py`). A carried C=C/C=N stereo bond was restored to DOUBLE
  unconditionally; when PuLP had relocated the double bond, that made an endpoint
  over-valent and every downstream `SanitizeMol`/`MolToSmiles` raised, so generation
  produced nothing. The promotion is now applied only when it keeps the molecule
  valence-valid, degrading to the documented "leave it and skip the constraint"
  behavior otherwise. Recovers `FIXYER`, `EDOFUB`, `EDOGEM`, `ZIHGEE` (all carry a
  `/C=C/` or `/C=N/` whose forced promotion over-valenced a carbon) to clean
  round-trips; `PILWUC` now generates too but reclassifies to `string_mismatch`.

* **A generation-internal wall-clock budget** (`embed_time_budget`, wired to the
  existing per-molecule `timeout`). The FF-only attempt loop had no time bound --
  `timeout` was consumed only by the ASE optimizer -- so a molecule whose embed never
  validated ran the full `max_attempts` (250) budget: `ZIHGEE_comp_0` took ~1696 s to
  return nothing. The loop now stops at the budget (checked between attempts, so a
  molecule that *does* embed is never interrupted) and returns whatever it has,
  turning a pathological non-terminating case into a fast, honest failure. The bound
  is `budget + one in-flight attempt`; a single pathological *attempt* is still
  covered by the harness `--mol-timeout` SIGKILL.
