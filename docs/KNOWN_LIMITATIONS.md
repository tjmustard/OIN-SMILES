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

## sp3 / heteroatom atom-stereo (`atom_stereo`)

The forward encode (`XYZToSMILES.convert`) and the generator's contract-mol re-encode
(`build_contract_mol` -> `get_oin_string`) perceive sp3/heteroatom chirality independently,
and the `atom_stereo` registry rows are where the two `@`/`@@` sets disagree. The combined
v0.3.7 R5 fix -- clearing stereo the encode never specified, plus re-orienting specified
centres on the metal-free fragment with `rdCIPLabeler` -- takes this class from **0/25 to
18/25 full round-trip successes**, and there are **no `@`/`@@`-disagreement residuals left**:
the remaining 7 rows all have their `@`/`@@` resolved and now fail on a *different* class
(donor-H atom count, or geometry RMSD / no conformer), i.e. they leave `atom_stereo`.

### Fixed: spurious stereo the forward encode never specified

- **Invented sp3 tags (e.g. `KAPCEM`, 0 `@` in -> 4 out).** The OIN leaves a ligand's sp3
  centres unspecified -- the fully-sanitised `get_tmc_mol` perceives them as non-stereogenic in
  the metal-bound complex. `build_contract_mol` then ran `AssignStereochemistryFrom3D` on the
  single embed conformer, which stamps a tag on *every* chiral-looking centre, inventing an `@`
  the input never emitted. `build_contract_mol` now records which sp3 centres the **parsed OIN
  template** specified and clears geometry-derived tags on any non-N/P centre the OIN left
  unspecified. Specified centres (a superset of the `sp3_stereo_targets` perceive-then-flip
  carry set) are untouched, so a genuine diastereomer is never masked. Guard:
  `tests/unit/test_heteroatom_atom_chirality.py`.
- **Spurious `-SF5` octahedral stereo (`MEDHUB`).** A pentafluorosulfanyl sulfur is octahedral
  with five identical terminal fluorines, so it is achiral -- but `AssignAtomChiralTagsFromStructure`
  stamps a `CHI_OCTAHEDRAL` tag from geometry and neither legacy nor modern RDKit perception
  (`FindPotentialStereo` calls it `Specified`) reduces the equivalent F, so the forward encode
  emitted a spurious `[S@OH..]` the generator correctly drops. `CIPAssigner` now clears
  high-coordination (`CHI_OCTAHEDRAL`/`CHI_TRIGONALBIPYRAMIDAL`) tags on non-metal centres whose
  stereochemistry rests on a set of identical terminal ligands.

### Fixed: wrong handedness at a metal-/η-adjacent specified sp3 centre

A specified sp3 centre bonded to a metal-bound ligand (e.g. `AHEBEV`'s benzylic carbon on an
η6-arene, `DAXJUI`, `KEBBUO`'s seven spiro-siloxane Si) round-tripped with the **wrong
handedness**: its CIP label *flips* between the metal-present contract mol (`R`) and the
metal-free fragment the `@` is actually emitted from (`S`), so the generator's metal-present
`sp3_stereo_targets` flip loop -- which compared the metal-free template label against the
metal-**present** contract label -- mis-oriented it. `build_contract_mol` now stamps the
metal-free template's `rdCIPLabeler` label (`_OIN_CIPCode_SP3`), and `ChiralityRecoveryUtility.recover`
verifies-and-flips the centre against it on the metal-free fragment (mirroring the Zone-A P
lone-pair branch; no-op on the forward-encode path, which never stamps the property). The label
is taken **aromatic-preserving on a fresh re-parse** of the template SMILES
(`_template_sp3_label` / `_reparse_aromatic_cip_label`, `SANITIZE_ALL ^ SANITIZE_KEKULIZE` with
an atom-map probe): rdCIPLabeler gives opposite R/S for a carbon bonded to an aromatic haptic
ring (η-Cp, and *fused* indenyl/fluorenyl -- `BABWAD`, `KAGXUM`, `NOSGAD`) depending on whether
the ring is aromatic vs kekulized, and the processed template object carries a corrupted
aromatic state, so the label is taken on a clean re-parse in the aromatic convention the
emitted fragment uses. **Both** the stamp and recover()'s comparison read the label this way --
using the re-parse on only one side flips a P-and-arene-adjacent centre (`BEPXEA`). The same
re-parse is applied to the **Zone-A P donor** lone-pair label (`_template_lp_label` +
recover()'s lone-pair branch), fixing a diphosphine donor bonded to an aromatic arm (`GUXPIA`,
whose `@` now matches -- it fails instead on a separate donor-H atom count). This carries the
18 successes above -- including cases previously mis-read as "wrong diastereomer," which the
RMSD gate confirms are geometrically correct (`DAXJUI` rmsd 0.52).

### Documented residuals

- **Stereo resolved, blocked by another class (leaves `atom_stereo`):** `KAPCEM`, `XENNIO`,
  `GUXPIA` (donor-H atom-count, S1 domain -- the `@`/`@@` now matches); `EJUKUQ`, `FADSAE` (High RMSD --
  the generated geometry is genuinely distorted, a conformer-quality/R2 issue); `GAKZOK`,
  `QOFTOU` (no conformer at all -- see the `no_conformers` entry; `QOFTOU` also builds rac/meso
  non-deterministically).

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

---

## Curated metal–ligand bond lengths cover a validated metal subset

The default generator places each σ donor at its FF-clean scan target
(`generator3d/clean_geometry.py::_binding_distance`). Historically this was a generic
Cordero covalent-radius sum, which systematically *over*estimates dative metal–ligand
bonds. A hand-curated per-`(metal, ligand)` table
(`generator3d/bond_lengths.py::BOND_LENGTHS`) now supplies realistic σ-donor
distances — but it is **applied only to `ENABLED_METALS`**, a dataset-validated subset
(`Ni, Pd, Pt, Zn, Cd, Hg, Ag`), not to all 18 metals the table carries.

- **Why a subset.** Validated against real input geometries (coordination-sphere mean
  RMSD, ~6 molecules/metal × 3 seeds), the table strictly improves the late/post-TM
  d⁸/d¹⁰ metals (median RMSD −0.015 to −0.084 Å) but *regresses* several
  early/mid-transition-metal buckets. The regressions are chemical, not stochastic:
  some entries encode a shorter bond mode than the dative bond present — `Ti–O 1.80`
  and `V–O 1.60` are metal-**oxo** (M=O) distances, and the early-TM `M–C` values are
  shorter than a real σ `M–C`. The generator cannot tell an oxo from an alkoxide, or a
  carbene from an alkyl, at this seam, so those metals keep the covalent sum. A few
  metals with a strong-but-inconsistent median (`Rh`, `Ir`) or a flat one (`Cu`, `Re`,
  `Fe`) are also left off. For every non-enabled metal the output is **byte-identical**
  to the pre-table generator — the change is strictly additive. Owned by
  `generator3d/bond_lengths.py`; expand `ENABLED_METALS` only with the same per-metal
  RMSD validation, never by hand.
- **The table is duplicated.** `generator3d/bond_lengths.py::BOND_LENGTHS` is a
  verbatim copy of `generation/molassembler_adapter.py::_BOND_LENGTHS` (the legacy
  backend that owns the original). It is copied, not imported, because importing that
  module pulls in Molassembler and the whole legacy generation backend. The copy is
  drift-guarded — `tests/unit/test_bond_lengths.py` asserts the two stay byte-equal, so
  an edit to either that is not mirrored fails CI (TD-005: a hand-copied constant is how
  Sc/Y once went missing).
