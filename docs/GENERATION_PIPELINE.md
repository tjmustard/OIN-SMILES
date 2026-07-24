# OIN-SMILES → 3D Structure: the Generation Pipeline

How the generator turns a 1D OIN-SMILES string into a 3D XYZ structure, step by step,
and how the result is validated. This documents the **default** path on current `main`;
where a step arrived in a specific wave it is marked (MetalloGen = the pre-existing vendored
engine; A5 = v0.4.3; SL0/SL1/SL4 = v0.4.4).

**One-line model:** parse the OIN string → build the metal-complex object *directly from the
parsed OIN* → distance-geometry-embed the whole complex against a dummy metal, clean it with a
constrained force field, and **stop at the first conformer that provably re-encodes to the exact
input** (falling back to geometry-based selection if none does) → rebuild the bonded mol and
emit XYZ. Correctness is proven by an independent XYZ→OIN re-encode compared with the fac/mer
canonical key.

---

## 0 · Configuration — what "default" resolves to

- **Entry point:** `OIN3DGenerator.generate(oin_string)` (`generation/engine.py`), wrapped by
  the public `SMILESToXYZ.convert`.
- **Resolved default knobs:** `timeout=300`s, `ensemble_size=10`, `optimizer="ff"` (no
  g-xTB/MACE relax), `seed=42` (fully deterministic).
- **Active-by-default behaviours:**
  - **OIN-direct assembly** — build the complex straight from the parsed OIN, no m-SMILES
    bridge (`OIN_DIRECT_DG` on; v0.4.4).
  - **vdW-clash acceptance gate** — reject sterically clashing conformers (`OIN_VDW_ACCEPTANCE`
    on; A5).
  - **Accept-first early-exit** — stop at the first conformer that round-trips (`OIN_EARLY_EXIT`
    on; SL1).
- **Opt-in-only levers (OFF by default):** rigid Kabsch placement / deterministic winding
  construction (`oin_direct`, SL2), greedy ordered placement (`greedy`, SL3), stretched-bond
  acceptance metric (`stretched_bond`, SL1), g-xTB/MACE optimizer.

---

## 1 · Parse the OIN string → `ParsedOIN`  *(MetalloGen-era)*

- `generation/oin_parser.OINParser.parse` reads the inline V3.6 string
  (e.g. `[Pt_SPL].N{0}.[Cl]{1}`).
- Extracts:
  - the **fragment SMILES list** split on `.` — the metal is always `fragments[0]` (a
    load-bearing canonical invariant);
  - the **geometry code** from the metal bracket (`_SPL`, `_OCT`, …);
  - a list of **`OINVector`** objects — one per binding atom — carrying fragment index,
    in-fragment atom index, **coordination slot number** `{n}`, and **winding** marker
    (`>`/`<`) for haptic rings.
- Each slot's ideal direction is looked up in a **generation-side template table**
  (`TEMPLATES`), kept deliberately separate from the encoder's table.

---

## 2 · Build the metal-complex object — directly from the OIN  *(default: v0.4.4; fallback: MetalloGen-era)*

The default path builds the internal `MetalComplex` **straight from `ParsedOIN`**, keeping the
metal + all ligand information (slots, winding, metal chirality) in one representation with no
lossy string round-trip. The old `metal|lig|…|geo` m-SMILES string is retained **only as a
fallback**.

- **Per-fragment preparation** — `metallogen_adapter._prepare_ligand_fragments` produces one
  `(mapped_smiles, winding)` spec per ligand (ordered by ascending first coordination slot).
  For each fragment it does the same chemistry the m-SMILES path always did:
  - **Slot → coordinate-vector matching (isomerism-preserving):** each binding atom is tagged
    with the atom-map number of the **nearest** MetalloGen coordinate slot vector. This
    nearest-vector match is what keeps cisplatin *cis*, and fac vs mer distinct.
  - **Donor hydrogen reconciliation:** a bare `C`/`N`/`O`/`S` donor sheds a phantom implicit H
    (anionic/dative donor); a bracketed H-bearing donor keeps its H.
  - **Kekulisation fix-ups** for neutral-radical aromatic rings (Cp) that can't kekulise
    normally.
  - **Winding capture:** for an eta ring whose heading atom carries a `>`/`<` marker, the
    canonical ring order + star atom + winding char are recorded, so the target winding can be
    reconstructed after RDKit canonicalisation.
- **Direct assembly** — `om.get_om_from_parsed` builds the `MetalComplex` from those specs:
  per-ligand connectivity (recovering C=C stereo, sp3 chirality, adding explicit Hs, and
  `binding_infos = [donor_atom_indices, slot]`), the geometry, and spin multiplicity (a
  weak-field high-spin heuristic, inert on the FF path). For a **load-bearing** haptic ring it
  attaches the winding target to the ligand (`Ligand.winding`).
- **What the carried extras do (and don't) in the default path:**
  - `Ligand.winding` is **carried but unused** by the default DG embed — only the opt-in rigid
    placer (`_place_haptic`) reads it. So haptic winding is still recovered by the *search*
    (§6), exactly as before. (Constructing winding from it was tried and regressed — see
    `docs/DIRECT_DG_VALIDATION.md`.)
  - **Metal `@SPn` chirality** rides along in the representation but is currently **inert for
    round-trip**, because the encoder cannot yet emit a reproducible metal stereo descriptor.
    It is groundwork for future encoder work.
- **Why direct, not m-SMILES:** removing the parse → serialize → parse round-trip means
  OIN-format changes reach 3D generation without a translation layer. An A/B matched the
  m-SMILES path byte-for-byte per molecule (0 regressions, 0 gains) — see
  `docs/DIRECT_DG_VALIDATION.md`.
- **Fallback** — if direct assembly raises on an edge case, `generate()` drops to
  `convert_parsed_to_msmiles` + `om.get_om_from_modified_smiles` (the winding-lossy m-SMILES
  path) rather than hard-fail. A fragment with **no** coordination slot (outer-sphere
  counterion / solvent) raises `UncoordinatedFragmentError` on either path.

---

## 3 · Build candidate 3D conformers — the embed loop  *(MetalloGen-era core; SL1/SL4 controls)*

`generator3d/__init__.generate_3d_structures` iterates `(scale, option)` combinations
(default `options=[0,1,2]`, `scales=[0.8…1.2]`) up to a budget, filling a pool — but early-exit
(§3e) usually cuts it short.

### 3a · Per-attempt: dummy-metal distance-geometry embed  *(MetalloGen-era)*
- `embed.get_embedding` swaps the real metal for a **dummy main-group centre** and places a
  **dummy atom** on each binding site.
- A single **RDKit distance-geometry embed** solves the *whole complex at once*, pinned by a
  **CoordMap**: metal → origin, each donor → `direction_vector[slot] × bond_length`.
  - Donors are **pinned** to their slot coordinates; ligand **bodies** are grown from a
    seeded-random start by the distance-geometry solve (satisfying bond/angle constraints).
  - `useRandomCoords=True` is the primary start; a metric-matrix retry rescues `-1` failures.
  - carried **E/Z and sp3 chirality** are enforced as distance-geometry constraints.
  - haptic binding triggers a `scales_for_haptic=[0.4…0.7]` sweep (the extra per-attempt cost
    that makes eta slow).

### 3b · Bond-length repair + validity gates — `_finalize_positions`  *(MetalloGen-era + A5)*
- **Bond-length patch:** over-long bonds from the atom substitution are pulled back to physical
  length.
- **Acceptance criteria (a conformer must pass all):**
  - **min interatomic distance** — no two atoms closer than a hard floor;
  - **inter-ligand collapse ratio** — no ligand pair fused below a covalent ratio;
  - **undesired-bond check** — perceived adjacency must match the input adjacency;
  - **vdW-clash count** *(A5, ON)* — reject any conformer with a non-bonded, non-geminal atom
    pair inside van-der-Waals contact (kept only as a scored fallback).

### 3c · Constrained force-field cleanup  *(MetalloGen-era)*
- `clean_geometry.TMCOptimizer.clean_geometry` runs a constrained MMFF/UFF minimisation with the
  metal–donor distances **pinned** (spring, not weld), relaxing ligand backbones without letting
  donors drift off their slots; curated metal–donor target lengths for enabled σ-donor pairs.

### 3d · Stereo filter  *(MetalloGen-era)*
- The cleaned conformer is checked for correct **E/Z** (dihedral) and **sp3 handedness** (signed
  volume) against carried targets; wrong-stereo conformers are diverted to `stereo_rejects` (a
  last-resort fallback only).

### 3e · Accept-first early-exit — **SL1, default-ON**  ← the key selection behaviour
- Each freshly-cleaned conformer is tested:
  - **independently re-encoded** — its XYZ is written and fed back through the *real* encoder
    (`XYZToSMILES().convert`), not the generator's own shortcut;
  - the re-encoded string's **fac/mer canonical key** (`canonical_roundtrip_key`, SL0) is
    compared to the input OIN's key; a cheap contract-mol pre-filter rejects obvious mismatches
    first.
- **On the first key match: stop building the pool and return that conformer** — skipping
  dedup, ranking, and the optimizer pass (§4).
- **Non-regressive by construction:** if the whole pool is built and nothing matches, it falls
  through to ordinary selection (§6).

### 3f · Loop control  *(MetalloGen-era + SL4)*
- **Wall-clock deadline** — stop between attempts when the per-molecule budget is exceeded.
- **No-progress cutoff** *(SL4)* — stop when N consecutive attempts add nothing (the gate is
  rejecting everything).
- **Structural-error surfacing** — if *every* attempt fails with the same structural error
  (e.g. an over-valent dative donor), raise a typed `StructuralAssemblyError`.

---

## 4 · Pool finishing — *only if early-exit did not fire*  *(MetalloGen-era + A5)*

- **Sort by UFF energy**, then **deduplicate** by heavy-atom RMSD + energy.
- **Re-rank by whole-complex vdW clash** *(A5)* — least-clashing first.
- **Optimizer pass** — with `optimizer="ff"` (default) this is a **no-op**; only `xtb`/MACE
  would relax the pool here.
- Return the top `num_conformers` (`max(ensemble_size, 5)`, or `16` when winding is requested).

---

## 5 · Coordination geometry / winding search  *(MetalloGen-era)*

- When the OIN encodes **eta-ring winding**, the pool is widened (`ETA_SELECT_POOL=16` + a
  doubled UFF pre-pool) so the requested ring face / diastereomer is actually sampled. Winding
  is thus recovered by **searching** the pool for a conformer that re-encodes with the right
  winding — this is why eta molecules build a wider pool.

---

## 6 · Select the returned conformer — `_select_by_geometry`  *(MetalloGen-era + SL0/SL1)*

- **Accept-first pass first** *(SL1)*: scan the returned conformers for one whose independent
  re-encode matches the input's fac/mer key; if found, return it (belt-and-suspenders with §3e).
- **Otherwise, the geometry ranking:**
  - keep only conformers whose coordination sphere **classifies** as the target geometry code;
  - if winding is requested, prefer the conformer whose re-encoded winding multiset matches;
  - rank survivors by `(vdW clash, template-fit RMSD, energy)`.
- Returns the chosen conformer plus its **contract mol**.

---

## 7 · Build the contract mol + emit XYZ  *(MetalloGen-era)*

- `build_contract_mol` reconstructs an RDKit molecule from MetalloGen's connectivity + coords:
  - **bond orders, aromaticity, formal charges, and stereo** transferred per-fragment from the
    OIN template SMILES via a **slot-anchored substructure match** (robust to non-sanitisable
    fragments like C#O or charge-less Cp);
  - metal→donor bonds set to **DATIVE**; chelate-locked E/Z stripped to match the encoder.
- `get_xyz_string` emits the final XYZ block. Result: `GeneratedStructure(xyz, mol)` — `.xyz` is
  the deliverable; `.mol` is the bonded contract mol (or `None` on an eta fallback).

---

## 8 · Validation

### 8a · In-generation acceptance (does a conformer qualify?)  *(MetalloGen-era + A5 + SL1)*
- Every conformer must pass, in order: **min-distance**, **inter-ligand collapse**,
  **undesired-bond adjacency**, and the **vdW-clash gate** (§3b).
- It must reproduce the carried **E/Z and sp3 stereo** (§3d).
- Under the default early-exit, it must additionally **re-encode to the exact input fac/mer
  key** to be accepted-and-returned (§3e) — generation validates against the round-trip
  objective itself.

### 8b · Round-trip validation (the external proof of losslessness)  *(MetalloGen harness + SL0 + SL4)*
- The generated XYZ is re-encoded **independently**, XYZ→OIN, via the full encoder
  (`get_tmc_mol` → `CIPAssigner` → aligner → serialiser).
- **String equivalence** is judged with the fac/mer-aware canonical key
  (`canonical_roundtrip_key`, SL0):
  - folds away benign relabelling (a rotation of the same isomer),
  - but distinguishes genuine isomers (fac vs mer) — the honest pass/fail signal.
- **Coordination-sphere RMSD** is a **diagnostic** *(demoted from a hard gate by SL4)* — a
  geometric sanity number, no longer able to fail a string-correct structure.
- **Bucketing** (`tools/roundtrip_bucket_report.py`, SL0) classifies each result:
  `byte_exact` / `key_equal` (benign) / `facmer_divergent` (real isomer error) / `structural` /
  `hard_fail` / `encode_fail` — the metric all accuracy work is measured against.

---

## Appendix · Provenance legend

- **MetalloGen-era** — the vendored dummy-metal + CoordMap DG embed + constrained-FF cleanup +
  contract-mol pipeline (v0.4.0–v0.4.2 base).
- **A5 (v0.4.3)** — vdW-clash acceptance gate default-on; engine `optimizer="ff"` default.
- **SL0 (v0.4.4)** — fac/mer-aware canonical key + bucket report (the measurement instrument).
- **SL1 (v0.4.4)** — accept-first early-exit (default-on); stretched-bond metric (opt-in).
- **SL4 (v0.4.4)** — RMSD demoted to a diagnostic; no-progress cutoff.
- **v0.4.4 direct-assembly default** — OIN-direct build replaces the m-SMILES bridge as the
  default (`docs/DIRECT_DG_VALIDATION.md`); m-SMILES retained as fallback.
- **Opt-in (not in the default path):** rigid Kabsch winding construction (SL2 `oin_direct`),
  greedy ordered placement (SL3), g-xTB/MACE optimizer.
