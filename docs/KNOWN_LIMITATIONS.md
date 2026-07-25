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

### Fixed in v0.4.2 (S6b): spurious donor-S, sulfonimidoyl S, quaternary N+, alkene/carbene-donor C

The full-dataset sweep re-surfaced two `@`-disagreement subsets on the v0.4.1 (`c7edeeb6`) baseline;
S6b closes them on the shared encode path (`ChiralityRecoveryUtility.recover`, which both
`XYZToSMILES.convert` and the contract-mol re-encode funnel through, so a fix is symmetric):

- **Spurious high-coordination donor-S (`[S@SP3]`/`[S@SP1]`/`[S@TB9H]`; `BAZMOH`, `HUGSEI`, `LUSKIV`,
  `YUMPIH`, `CIDDAU`).** A metal-donor thioether/ring S gets a permutation chiral tag from
  `build_contract_mol`'s `AssignStereochemistryFrom3D` (metal-present geometry) that survives
  fragmentation because legacy `AssignStereochemistry(cleanIt=True)` does not scrub a pre-set
  permutation tag -- the input crystal geometry never produced it. `recover()` now clears a chiral
  tag on an S that is **not** a genuine stereocentre on the metal-free fragment
  (`FindMolChiralCenters(useLegacyImplementation=False)` -- a divalent thioether S is absent; a genuine
  sulfonimidoyl S(VI) is present and kept). `CIDDAU`'s `@` now matches; it fails instead on a separate
  `[SH]` donor-H atom count (S1 domain).
- **Genuine sulfonimidoyl S(VI) inverted (`JEKQAS`, `REPZUJ`, `ZORCOA`).** The centre's metal-donor
  O becomes a radical `[O]` in the fragment, which `rdCIPLabeler` refuses to rank, so
  `_reparse_aromatic_cip_label` returned no label and the `_OIN_CIPCode_SP3` re-orientation was
  silently skipped. The reparse now **fills a metal-stripped donor's open valence with H** (skipping
  aromatic atoms) so the CIP is computable, in the same convention `_template_sp3_label` reads.
- **Genuine quaternary ammonium N+ dropped (`POYJIX`, `[N@@+]` -> `[N+]`).** `build_contract_mol` did
  not stamp N, so `recover()`'s 4-neighbour no-`_OIN_CIPCode` fallback cleared it. A **degree-4** N+ is
  now routed through the same metal-free `_OIN_CIPCode_SP3` re-orientation as C/Si/S (a trivalent
  amine N stays unstamped -- see below).
- **Genuine sp3 C bonded to a metal-bound alkene/carbene donor inverted (`ORIHUU`, `XILZID`).** The
  donor carbon is valence-deficient in the fragment, so its template-vs-fragment CIP diverged and the
  re-orientation mis-fired. `_template_sp3_label` and `_reparse_aromatic_cip_label` now read the label
  through the **same fill-first reparse**, so the donor is normalised identically on both sides.

### Documented residuals

- **Trivalent-N inversion is unrepresentable (`JUCCUH`, `[N@@H]` -> `[NH]`).** RDKit clears a
  trivalent amine `[N@]` as a non-stereogenic nitrogen inversion, so the generated re-encode cannot
  carry it. Zone-A / backbone trivalent N stereo remains **deferred** (needs an out-of-band marker);
  this is the same limitation R5 documented for Zone-A N.
- **Macrocyclic multi-Zone-A-P relative configuration (`WEDYOU`).** A 1,4,7-triphosphacyclononane on
  Fe whose two stereogenic P donors round-trip with a swapped **relative** configuration
  (`[P@H]…[P@@H]` vs `[P@@H]…[P@H]`; canonical keys genuinely differ). The embed builds one relative
  diastereomer non-deterministically -- a generator stereo defect analogous to `QOFTOU`'s rac/meso,
  not a per-centre carry gap the lone-pair re-orientation can fix; **deferred**.
- **Stereo resolved, blocked by another class (leaves `atom_stereo`):** `KAPCEM`, `XENNIO`,
  `GUXPIA` (donor-H atom-count, S1 domain -- the `@`/`@@` now matches); `EJUKUQ`, `FADSAE` (High RMSD --
  the generated geometry is genuinely distorted, a conformer-quality/R2 issue); `GAKZOK`,
  `QOFTOU` (no conformer at all -- see the `no_conformers` entry; `QOFTOU` also builds rac/meso
  non-deterministically); `CIDDAU` (`[SH]` donor-H atom count, S1 domain -- the `@` now matches).

---

## Borane and carborane clusters (`carborane_unsupported`)

`get_lig_mol` (`utils/xyz2mol.py:426`, calling `AC2mol` at `:451-462`) cannot perceive
bond orders for electron-deficient polyhedral cages (`OSENOR`, and the
`carborane_unsupported` bucket). It fails in the **forward** encode with a specific
message (`get_lig_mol failed for ligand fragment #... (SMILES: '[H]B1[B-]2...')`), and
the charge sweep does not rescue it: the 3-center-2-electron bonding these cages need is
outside the **two-centre** bond-order model `AC2mol` implements.

**Layer that must change:** the OIN notation itself — a cluster convention (an eta-like
multi-atom unit, or a pseudo-atom for the cage) — not the encoder. This is deferred until
a design exists; it is a notation-design gap, not a bug.

**Members (snapshot 2026-07-14):** 36 on the `c7edeeb6` floor, 92 in the current
(mixed-provenance, growing) backlog. Full list and root-cause detail:
`spec/handoffs/v0.4.2/wontfix-carboranes.md`. Regenerate with
`tools/classify_failures.py` on a **copy** of `results-v0.4.0/` (class
`carborane_unsupported`).

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
gap for that ligand set. (v0.4.2 note: adding the CN-9 `TCT` geometry code made
`XERTUK_comp_3` *encodable* as `[Y_TCT]` -- previously `g:NON` -- but its 104-atom
ligand still will not embed; a geometry code being supported does not guarantee the
builder can realize a large, crowded ligand set.)

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
- **Re-validated in v0.4.2 (S7): the subset stands.** A paired median-of-deltas coordination-sphere
  RMSD A/B (≥10 seeds, **full** conformer pool — collapsing the pool fabricates results) with a Zn
  positive control (−0.025 Å, reproducing the landed win, so the harness is trusted) confirmed no
  candidate metal earns inclusion: Ru comes back **+0.009 Å** (no benefit — the swept scale factor
  already compensates the covalent overestimate) and W/Mn/Co are flat. `bond_lengths.py` ships
  **unchanged**; the `high_rmsd` bucket is the FF-only geometric floor, not a bond-length bug (below).
- **Provenance of the table.** `generator3d/bond_lengths.py::BOND_LENGTHS` is the single
  authoritative copy. It originated as a verbatim copy of the legacy Molassembler
  backend's `_BOND_LENGTHS` (that backend was removed in v0.3.7). The values are
  drift-guarded — `tests/unit/test_bond_lengths.py` asserts them against the frozen
  legacy table, so an unmirrored hand edit fails CI (TD-005: a hand-copied constant is
  how Sc/Y once went missing).

---

## Donor hydrogen count -- the nitride/ammine notation ambiguity (`donor_H_atom_count`)

An OIN string carries no explicit hydrogen on a de-bracketed organic-subset donor, so the
generator re-materializes donor H with `Chem.AddHs` (`generator3d/ligand.py:63`) and the
adapter's H-reconcile decides which bare donors are 0-H X-type
(`generation/metallogen_adapter.py:164-222`). One case there is **genuinely undecidable at
the notation layer** (documented in-code at `:207-219`): a **bare, heavy==0 nitrogen** donor
-- a nitride `[N]` bound only to the metal and an ammine `[NH3]` -- can serialize to the same
`N{n}` token, and the reconcile lets **the ammine reading win** (it keeps the three H). A
nitride so written round-trips with +3 spurious H.

**Layer that must change:** the OIN format (`oin/inline.py:265-317`). The encoder already
disambiguates a 0-degree nitride to a bracketed `[N]{n}` (gated on `GetDegree()==0`,
`inline.py:300-306`); closing the residual needs that marker extended to every bare-N case
the generator would otherwise read as ammine, verified not to regress real ammine donors. It
is a notation-design change, not a generator bug.

**Members / current status (snapshot 2026-07-14).** `donor_H_atom_count` is 82 on the
`c7edeeb6` floor / 202 in the current backlog, and **S1 (`feature/roundtrip-donor-h`) owns the
fixable-vs-notation split** of that class -- the notation-limited subset it routes here is
finalized when S1 lands; this phase then absorbs it. Two honest caveats on the accounting as
it stands today, so a later reader is not misled:

- A fresh `classify_failures.py` run isolates **no** `donor_H_atom_count` row by the
  lone-`N{n}` nitride/ammine signature -- so this ambiguity is presently a **latent** format
  limitation with no clean dataset trigger, not a populated bucket.
- The `+3`-atom donor-H rows are **not** this class. They are dominated by
  eta-arene / `[CH]`-radical ring-H re-materialization -- e.g. `AJIJUY_comp_0` re-encodes
  **byte-identically** (`smiles_1 == smiles_2`) yet the generated 3D carries +3 H, and the
  fragment contains no nitrogen at all. Those belong to S1's generator-fixable / eta H-count
  work, not to the nitride/ammine notation limit.

---

## FF-floor high RMSD and `--quick` timeouts are harness artifacts, not accuracy failures

Two of the largest "failure" buckets in the sweep are **not** round-trip accuracy defects; they
are properties of the benchmark harness and the FF-only generation decision.

- **`high_rmsd` (36 floor / 112 current).** The 1.0 Angstrom RMSD gate
  (`tools/test_dataset_roundtrip.py:173`) runs **after** the OIN string round-trip has already
  matched -- so every `high_rmsd` row is a **chemically-correct** round-trip that failed only on
  geometric tightness. Under the wave's FF-only decision (no `xtb` binary -> the `g-xtb` optimizer
  warns and returns the FF geometry unchanged), coordinate-sphere RMSD is systematically inflated
  for any metal outside `generator3d/bond_lengths.py::ENABLED_METALS` (`{Ni,Pd,Pt,Zn,Cd,Hg,Ag}`),
  which falls back to a covalent-radius sum that **over**estimates dative bonds (see the curated
  bond-length section above). Do **not** loosen the gate to "fix" these; the real, bounded win is
  extending `ENABLED_METALS` with per-metal-validated bond lengths (S7's `feature/roundtrip-metrics`).
- **`timeout` (339 floor / 948 current).** A `--quick` **labeling** artifact: the accumulator runs
  `ff_params_fast = {uff_pool_size: 2, max_attempts: 10}` with a 30 s hard-kill
  (`test_dataset_roundtrip.py:616-617`, SIGKILL at `:253/:302`), so a slow molecule is filed
  `timeout` regardless of *why* it is slow. A v0.4.2 (S7) full-budget confirmatory sample found most
  are **real valence `no_conformers`** the 30 s kill masks -- they still fail at full budget, with the
  honest label -- not molecules a longer budget recovers. Either way they are not string-level
  accuracy defects; the harness now stamps `mol_timeout`/`rmsd_gate` and prints a failure breakdown so
  `timeout` / FF-floor rows are not counted as accuracy failures.

**Layer:** the harness thresholds / the FF-only decision -- not the encoder or the generator. S7
owns the honest split (genuine vs artifact) and the harness relabeling; docs records that these
buckets are, by construction, not accuracy regressions.

### Full-quality generation can time out where quick mode succeeds (observed v0.4.4)

The v0.4.0-vs-v0.4.4 regression sweep (`REGRESSION_REPORT.md`) surfaced a *different* timeout
regime from the `--quick` labeling artifact above. Of the 1,000 v0.4.0-success guard molecules,
**11 no longer round-tripped under full quality** -- and every one is a 300 s generation
`TimeoutException` (10 UFF, 1 g-xTB), **not** a wrong-answer/notation regression (zero regressions
landed in `structural`/`facmer_divergent`/`encode_fail`). All 11 are medium-large (69--157 atoms)
and each finished in **6--29 s under v0.4.0 quick mode** (`uff_pool_size=2, max_attempts=10`, 30 s
cap) but exceeded the 300 s cap under full quality (full UFF pool + direct-DG default) -- a
10--50x per-molecule slowdown (`TIVNAQ`/`XEXMEV` went 6--7 s -> >300 s). Members: `CALQEO`,
`DOSBUJ`, `DOZGIJ`, `HACDUM`, `LUZDEP`, `QOVSAV`, `QUSYUA`, `SOJMIQ`, `TIVNAQ`, `WIMBOL`, `XEXMEV`
(all `_comp_0`).

**Layer / status:** generator compute cost, not encoder/notation. Root cause is **not yet
isolated** -- full-pool size vs the direct-DG embed path -- so this is recorded as an observed
full-quality generation-time limitation, not a fixed defect. Open follow-up: re-run the 11 at
quick config (they each passed in <30 s at v0.4.0) to confirm no code defect, then profile the
full-pool + direct-DG path on these structures.

---

## Uncoordinated outer-sphere fragments (`UncoordinatedFragmentError`)

Outer-sphere counterions and uncoordinated solvent -- a free water, a borate/fluoroborate anion,
a perchlorate, a bare oxide -- are emitted as their own OIN fragments by the encoder but carry no
bond to the metal, so they have **no binding slot**. MetalloGen's `metal|lig1|...|geo` m-SMILES
has no way to express such a fragment, and generation raises `UncoordinatedFragmentError`
(`generation/metallogen_adapter.py:92`, raised at `:154`). Examples: `ATAGUZ_comp_0`
(`[O-]Cl` perchlorate + `[O-2]` oxide counterions on an `[Hg_LIN]`), `CAXZAD_comp_0` (a
perfluoro-tetraarylborate counter-anion).

**Layer:** representation. m-SMILES has no non-coordinating-fragment slot; round-tripping these
would need an OIN/generator convention for outer-sphere species. A representation limit, not a bug.

**Members (snapshot 2026-07-14):** at this snapshot **all** `gen_exception_other` rows are
`UncoordinatedFragmentError` -- 24 on the `c7edeeb6` floor, 72 in the current backlog. Floor set:
`ATAGUZ_comp_0, CAXZAD_comp_0, CEGBAU_comp_0, CUBDOT_comp_0, DEJHEF_comp_0, EJOPOJ_comp_0,
FOQBEV_comp_0, FUPVOC_comp_0, GEYRUA_comp_0, KOFLAV_comp_0, LOJGEW_comp_0, MAHTOE_comp_0,
NEBDAA_comp_0, OHAYIH_comp_0, SOLYEZ_comp_0, SOSJAM_comp_0, TIGDAO_comp_0, UHAMUM_comp_0,
VIDVUA_comp_0, XAXYEC_comp_0, XAXZIH_comp_0, XIFVAM_comp_0, XIQKOY_comp_0, YUMBEP_comp_0`.
Regenerate the full current set by grepping the report `error` for `UncoordinatedFragmentError`.

---

## Irreducible generator stereochemistry

Some stereochemistry the input carries cannot be reproduced deterministically by the FF-only
builder, so the round trip differs on stereo through no encoder fault.

- **Non-deterministic rac/meso construction (`QOFTOU`-class).** `QOFTOU_comp_0` builds a
  racemic/meso mixture non-deterministically at generation time -- the same seed can realize
  either diastereomer -- so its `@`/`@@` cannot be pinned. It already appears under the
  `atom_stereo` residuals above (its `@`/`@@` disagreement is resolved; it fails instead on
  no-conformer + this non-determinism).
- **Irreducible winding (`winding_flip` residual).** `winding_flip` (14 on the `c7edeeb6` floor /
  33 current) is **largely a `--quick` conformer-pool artifact**, not a builder defect: under the
  non-quick default generator (seed 42, pool 5) the correct coordinated face is realized for **12 of
  14** floor rows. v0.4.2 (S5) confirmed this and shipped **no winding code** (its only change was
  the CN-9 template); a small residual encodes a ring-winding sign the FF builder cannot reproduce
  deterministically. (For the fixed, geometry-inert cases -- Cp*/arene/BPh4- rings a proper rotation
  turns over -- see the eta-winding note in the project history; those are *not* residuals.)

**Layer:** generator geometry realization / builder stereo -- it cannot deterministically
reproduce input handedness at a metal-adjacent centre or for a winding whose sign a proper
rotation flips. Not an encoder or notation defect.

---

## Encoder injectivity blind spots (Y1 audit -- round-trip *false positives*)

Every limitation above is a round-trip **FALSE NEGATIVE** (the round trip FAILs, but the OIN
was fine -- usually the generator). This section is the opposite and more dangerous cell: axes
where two **genuinely distinct isomers encode to the same OIN**, so the round trip PASSES while
the notation is silently lossy. Found by the Y1 injectivity audit (`tools/injectivity/`), which
mirror-twins a structure and asks whether the encoder can still tell the enantiomers apart --
no 3D generator involved. Full write-up: `docs/INJECTIVITY_Y1_OVERVIEW.md`.

**Wave 2 update (feasibility):** all three axes are **recoverable from the input 3D** — the
encoder discards a signal it has, it does not lack one — so none is a *permanent* limitation
(`docs/INJECTIVITY_Y2_FEASIBILITY.md`). A per-axis descriptor recovered from the coordinates flips
between enantiomers on every fixture. P2 has an opt-in encoder emit; P1/P3 are deferred to the
v0.4.5 canonical-string work (they need a canonical, orientation-invariant ordering + generator
support before a token can be emitted safely).

- **Metal-centre Δ/Λ (@SP/@OH) -- KEY-BLIND (recoverable, deferred to v0.4.5).** `fac-Ir(ppy)3`
  and its enantiomer produce different *raw* strings, but only by non-reproducible slot renumbering,
  which the round-trip key deliberately folds (`oin/compare.py` `_METAL_STEREO_RE` +
  `_polyhedron_signature`). RDKit's `AssignStereochemistryFrom3D` *does* recover the Δ/Λ
  configuration (octahedral permutation 10 vs 8), so this is recoverable; it needs a canonical donor
  ordering (the v0.4.5 problem) + generator support to emit an `@OHn` token. Detail:
  `docs/INJECTIVITY_Y1_P1_METAL.md`.
- **Axial / atropisomeric chirality (biaryl, BINAP) -- ENCODER-BLIND by default; OPT-IN EMIT that
  now ROUND-TRIPS (Wave 2).** `R`-BINAP and `S`-BINAP encode to **byte-identical** OIN strings by
  default. The axis is recovered from the signed biaryl dihedral (`src/oinsmiles/oin/axial.py`);
  behind `OIN_EMIT_AXIAL` (default OFF) the encoder appends a *canonical* axial-sign token so the
  two diverge (`|ax:-|` vs `|ax:+|`), the round-trip key folds it (batch unaffected), and the
  generator honours it (axial-aware conformer selection + acceptance): measured **2/2 vs 1/2
  baseline**, so the axis survives the round trip. A stereogenicity gate keeps the encoder from
  claiming chirality for achiral symmetric-end biaryls. Still opt-in pending a broader A/B than
  one fixture pair. Guards: `tests/unit/test_axial_emit.py`,
  `tests/integration/test_axial_roundtrip.py`. Detail: `docs/INJECTIVITY_Y1_P2_AXIAL.md`.
- **Trivalent-N amine inversion on binding -- ENCODER-BLIND (recoverable, deferred to v0.4.5).** A
  *metal-bound* secondary amine is stereogenic only because the metal locks its 4th position; the
  encoder strips the metal, leaving a trivalent N that Zone-A clears (`core/chirality.py:722-727`; N
  is out of lone-pair-CIP scope, `:33-36`). Recoverable as the signed tetrahedral volume at the
  locked N (−9.4 vs +9.4); emitting needs a Zone-A carve-out + canonical ordering + generator.
  Distinct from the pendant-amine `[N@@H]→[NH]` residual under *Documented residuals* above, which
  is a generator loss. Detail: `docs/INJECTIVITY_Y1_P3_AMINE.md`.

**Layer:** encoder injectivity / notation completeness. Reproduce with
`PYTHONPATH=$PWD/src python -m tools.injectivity.report --probes`. Guards:
`tests/unit/test_injectivity_probes.py`, `tests/unit/test_config_oracle.py`,
`tests/unit/test_axial_emit.py`.
