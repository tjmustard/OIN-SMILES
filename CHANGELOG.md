# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.3.5] - 2026-07-08

### Added
- **Canonical symmetric-donor binding slot** (`utils/oin_aligner.py`, `utils/xyz2mol.py`): a monodentate ligand that binds through one of two resonance-equivalent atoms (e.g. a carboxylate's two oxygens, which differ as `=O`/`-O` in any single Kekulé structure) now always carries the `{slot}` marker on a canonically-chosen atom. Structures that differ only in which atom 3D bond perception happened to pick now encode identically, fixing a spurious round-trip "String mismatch". Guarded by `tests/unit/test_canonical_donor_binding.py`.
- **`tools/recalculate_oin_smiles.py`**: A utility to recalculate OIN SMILES strings from both the input XYZ and the generated XYZ structures for previously processed datasets. Updates the `summary_roundtrip.json` and `individual_reports` statuses if a codebase change causes a previously failed mismatch to now perfectly round-trip.
- **`tools/rebuild_summary.py`**: Rebuilds `summary_roundtrip.json` from the per-molecule `individual_reports` already on disk.
- **Dataset Roundtrip Tools enhancements**: Added `--quick`, `--continue`, `--rerun-failed`, `--random`, and `--mol-timeout` options to `tools/test_dataset_roundtrip.py` to allow robust, resumed, and time-bounded background processing of large datasets without hanging on pathological UFF geometries.

### Fixed
- **η³-allyl double-bond loss on round-trip** (`generation/metallogen_adapter.py`): `_flatten_template` built the connectivity-only query for `build_contract_mol`'s per-fragment substructure match but cleared only aromaticity and charge — not radical electrons. A ligand atom that binds the (stripped) metal is under-valent, so its template atom carried a radical, and `GetSubstructMatch` treats radical count as a match constraint — so bond-order/aromaticity transfer silently failed and the ligand was emitted all-single and de-aromatized (the η³-allyl double-bond loss, e.g. ABAZEK). Now clears radicals and normalizes H valence so the match succeeds and the allyl `=` is preserved by transfer; five dataset allyl cases (ABAZEK/ABETIK/ABETOQ/ACALOI/AGOVOK) now key-match. Guarded by `tests/unit/test_contract_mol_allyl_transfer.py`.
- **Zone-A chiral P donor stereo lost on round-trip** (`generation/metallogen_adapter.py`): a phosphorus binding the metal directly with a stereogenic lone pair encoded as `[P@]{0}` on XYZ→OIN but re-encoded achiral `P{0}` after OIN→3D→OIN (e.g. ACUWUT). `build_contract_mol` never populated the lone-pair path `recover()` needs for the donor, and the dative metal→P bond makes `AssignStereochemistryFrom3D` return `CHI_UNSPECIFIED`. The generated donor is now primed with `_OIN_CIPCode_LP` and a seeded chiral tag from `rdCIPLabeler` on the metal-free template (not the legacy `Chem.AssignStereochemistry` label, which disagrees for 3-coordinate P and would round-trip the wrong enantiomer), so `recover()`'s lone-pair verify-and-flip branch re-asserts the encoded handedness. Guarded by `tests/unit/test_zone_a_p_donor_stereo.py`.
- **`--quick` roundtrip crash from unsupported `max_attempts`** (`generator3d/__init__.py`): `--quick` mode passes `ff_params={"max_attempts": 10}`, which `generate_3d_structures` forwarded verbatim to `TMCOptimizer(**ff_params)` — but `TMCOptimizer.__init__` has no such parameter, so every molecule raised `TypeError: unexpected keyword argument 'max_attempts'` and crashed before any geometry was produced. `max_attempts` is now filtered out before constructing the optimizer while still capping the embedding retry loop.
- **UFF loop hang during 3D generation**: Fixed a bug in `generate_3d_structures` where an unrecognized UFF atom type (which immediately fails FF cleaning) caused the generator to blindly brute-force 250 random embeddings before giving up. Added support for a `max_attempts` override in `ff_params` (10 under `--quick`) that caps the embedding loop.

### Changed
- **g-xTB optimizer renamed `xtb` → `g-xtb`** across the user-facing surface: the `oin-smiles oin2xyz --optimizer` default, the `ASEOptimizer` method (which still accepts both spellings), and the dataset roundtrip harness tiers. The dataset harness now also short-circuits hard generation/verification failures and timeouts instead of escalating them to the slow g-xTB pass.

## [0.3.4] - 2026-07-07

### Fixed
- **Default g-xTB optimizer crash.** `ASEOptimizer.optimize()` unconditionally called `atoms.get_potential_energy()` after refinement, but the g-xTB path (a subprocess, not an ASE calculator) never attaches one — so every *successful* g-xTB optimization raised `Atoms object has no calculator`. The energy is now read from the ASE calculator only on the MACE path; the g-xTB path uses the energy parsed from `xtb` stdout. Guarded by a new regression test. (`generator3d/ml_optimizer.py`)
- **Non-functional charge/bond-order code paths.** `process.get_chg_and_bo` / `get_bo_matrix_from_adj_matrix` / `get_chg_list_from_bo_matrix` referenced a `frag`/`compute_scipy` module not vendored in this subset (would `NameError` if reached); they now route through the active PuLP solver (`compute_chg_and_bo_pulp`). Also fixed a missing `self` parameter on `Molecule.get_screening_result`, a `print_x` typo in the `Ligand` debug print, and a duplicate `get_electron_count` definition.

### Changed
- **MACE / PyTorch are now an opt-in extra.** `mace-torch` and the pinned CUDA-11.8 `torch` moved from hard dependencies to `pip install/uv sync --extra mace`. The default install is lightweight (FF + g-xTB, no `torch`), matching the documented "fast FF-only path". MACE optimizers (`mace-omol-*`) require the extra and still fail loudly if unavailable.
- Removed dead, unreachable, partially-ported code from the vendored MetalloGen engine (fragment-based charge/BO heuristics, `get_kth_neighbor_atom_list`, `scatter_molecules`, and the `frag`-based `Detect_EZ`/`Detect_RS`/`Detect_stereocenter` stereo path — OIN-SMILES handles stereo upstream).

### Added
- **Configurable g-xTB timeout with FF fallback.** A `timeout` (default 300s) threads through `OIN3DGenerator` -> `MetalloGenAdapter` -> `generate_3d_structures` -> `ASEOptimizer`; when `xtb` exceeds it, `subprocess.TimeoutExpired` is caught and generation falls back to the FF geometry instead of hanging (the old fixed 60s cap was too short for large complexes).
- **`tools/install_mace_weights.sh`**: idempotent downloader for the public `MACE-omol-0-extra-large-1024` checkpoint (ACEsuit GitHub release) that also registers `MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH` in `.env`. `models/mace/README.md` reworked accordingly (the extra-large model is freely downloadable; OMol25 remains Hugging-Face-gated).
- **`tests/unit/test_generator3d_units.py`** (32 tests): unit coverage for the vendored `generator3d` engine — `chem.Atom`/`Molecule`, `process` helpers, the frag→PuLP reroute, package helpers, and `ASEOptimizer` (incl. the g-xTB regression test).
- **Lint tooling & CI**: `ruff` added as a dev dependency (`[dependency-groups]`); the entire repo is now `ruff check`/`ruff format` clean (was 1817 findings). Added `.github/workflows/ci.yml` running lint + the unit and encoder suites on push/PR.
- **Packaging metadata**: `pyproject.toml` gains a proper description, `authors`, `license`, `keywords`, `classifiers`, and `[project.urls]`.
- **Docs**: `docs/OPTIMIZERS.md` (FF vs g-xTB vs MACE selection/install) and a `docs/README.md` index; README installation steps corrected (numbering + light-vs-MACE install).

### Removed
- Stray repository clutter: `os` (empty), `main.py` (hello-world stub), `REPORT.md`, `OpenSourceTMCBuilderReport.md`, a stray `verify_xyz_to_oin.py.fragment`, and four unrelated HACF-framework docs under `docs/`. The root HACF installer was renamed `install.sh` → `HACF-install.sh` to distinguish it from project installation (`uv sync`).

## [0.3.3] - 2026-07-06

### Changed
- **Default generation engine → MetalloGen; default optimizer → g-xTB.** `OIN3DGenerator.__init__` now defaults to `engine="metallogen"` and `optimizer="xtb"` (was `engine="legacy"`, `optimizer=None`), and `cli.py oin2xyz` rides that default. The MetalloGen backend is the better-validated generator (full MACE round-trip 25/25; eta-winding rac/meso, BDNN square-plane, and TiCat eta TET/TPY all fixed). g-xTB is fast and, unlike MACE, **degrades gracefully to FF** when the `xtb` binary is not on `PATH` (so the default never hard-fails). MACE (`mace-omol-0-extra-large-1024` / `mace-omol25`, higher accuracy, requires `mace-torch` + weights and fails loudly if absent) and the legacy SCINE Molassembler backend (via `engine="legacy"`, still the reference for Zone-A P stereo enforcement) remain opt-in.
- **`verify_roundtrip.py --optimizer` default → `xtb`** (was MACE) for fast iteration; pass `mace-omol-0-extra-large-1024` for the accurate sign-off.

### Added
- **g-xTB optimizer** (`generator3d/ml_optimizer.py`): a subprocess wrapper around the `xtb` binary (`xtb <struc> --gxtb --opt`) selected via `optimizer="xtb"`. Fast semi-empirical refinement of the FF pool; warns and falls back to FF if the binary is missing. Helper `tools/install_gxtb.sh` installs the binary; `tests/integration/run_optimization_grid.py` benchmarks FF/g-xTB/MACE across the fixtures.
- **`oin-smiles oin2xyz --engine {metallogen,legacy}` and `--optimizer`** (default `xtb`; accepts `ff`/`none`/`mace-omol-*`), giving a fast, no-torch, or legacy path from the CLI.

### Fixed
- **Legacy-specific real-generation unit tests pinned to `engine="legacy"`** (`test_winding_inertness.py`, `test_zone_a_p_genenforce.py`, `test_stereo_roundtrip_diagnostics.py`) so the default flip keeps the fast unit suite deterministic and heavy-optimizer-free while those Molassembler-only behaviors stay under test.

## [0.3.2] - 2026-07-06

### Added
- **Eta-ligand winding (rac/meso) round-trip fidelity**: the encoder now emits a winding marker per haptic slot — not just the first ring — measured against each ring's *actual* metal→centroid axis (`oin_aligner.py` `_permute_and_serialize` / `_determine_winding`, was the idealized template slot axis that flipped the 2nd ring under a distorted ansa bite). The MetalloGen generator honors it via winding-**multiset** conformer selection over a widened eta pool (`metallogen_adapter.py` `_eta_winding_multiset` / `_reencode_oin_fast` / `ETA_SELECT_POOL`), fixing the TiCat3/TiCat4 rac↔meso diastereomer swap (both now round-trip to the correct isomer). Generalized to N eta ligands and variable hapticity (η³ definite / η² degenerate); `verify_roundtrip.py` compares via a winding-canonical key (winding-stripped string + sorted multiset). Adds `TiCat5`/`TiCat6` fixtures and `test_eta_winding_generalization.py` (8 tests).
- **`oin_aligner.py`**: Added `classify_coordination_geometry()` (best-matching OIN geo code for a set of metal-centred donor vectors) and `coordination_geometry_fit()` (RMSD of the best donor-to-template assignment — the fit quality the classifier itself discards). Both wrap the existing discrete-geometry matcher.
- **`MetalloGenAdapter`**: Added `_select_by_geometry()` geometry-code-aware conformer selection. From the energy-ranked pool it keeps only conformers whose coordination sphere classifies as the requested geometry, then returns the tightest template fit (energy breaks ties). Haptic/η donors are gated out (donor count ≠ coordination number), making selection a deliberate no-op there and strictly non-regressive versus lowest-energy.
- **`generate_3d_structures`**: Added conformer deduplication over the FF pool via new `uff_pool_size`, `rmsd_threshold`, and `energy_threshold` parameters, plus a `calculate_heavy_atom_rmsd()` helper. `MetalloGenAdapter` surfaces these through `ff_params`.
- **`rmsd_utils.py`**: Added `_compute_robust_rmsd()` (anchor-pair candidate rotations → Hungarian assignment → Kabsch refine → ICP polish, floored by the greedy estimate) for the >5-atoms-per-element branch, fixing bent ansa-metallocene mis-pairing.
- **`test_geometry_selection.py`**: Added 18 unit tests covering the classifier, template-fit ranking, coordination perception, the haptic gate, and lowest-energy fallbacks.
- **`run_verification.sh`**: Added a `--limit N` pass-through to `verify_xyz_to_oin.py` and `verify_roundtrip.py`.
- **`tools/test_dataset_roundtrip.py`, `tools/test_uff_pool_size.py`**: Added dataset round-trip and UFF-pool-size sweep scripts.

### Fixed
- **`MetalloGenAdapter`**: Fixed the stochastic PdCl2-RR-BDNN failure where the generated Pd distorted from square-planar (`SPL`) toward trigonal-pyramidal, giving RMSD ~1.4 and a geo-code mismatch. Geometry-fit-ranked selection now returns the cleanest square-plane from the pool (BDNN: 5/5 round-trip PASS, RMSD 0.10–0.20; was intermittent). The prior `--ensemble-size` lever was a no-op under an optimizer — the pool is fixed at `pool_size` and only the lowest-energy conformer was used.
- **`generator3d/__init__.py`**: Fixed a `TypeError: '<' not supported between instances of 'NoneType'` crash that intermittently aborted generation when pool energies were unset, by making the energy sort None-safe and computing a final FF energy for each conformer.
- **`rmsd_utils.py`**: Fixed a `997` coordination-sphere false positive on ansa-metallocenes (TiCat1–4). The non-bonded ansa-bridge Si fell inside the distance cutoff and broke element-set equality; `calculate_tmc_rmsd` now drops from the distance-based input sphere any element absent from the bond-based generated sphere (bond sphere is donor ground truth).
- **`MetalloGenAdapter`**: Fixed TiCat2 Cp radical aromaticity loss by calling `Chem.RemoveHs(t, sanitize=False)` on sanitize-failed fragments, so the re-encoded OIN keeps aromatic `c1[cH]…` instead of kekulized `C1[CH]=…`.
- **`MetalloGenAdapter`**: Closed the stochastic TiCat1/3 `[Ti_TET]`↔`[Ti_TPY]` round-trip string drift by extending geometry-fit-ranked selection to haptic ligands. `_coordination_vectors` now reduces hapticity to centroid donors (new `_reduce_haptic_positions`: <1.6 Å transitive clustering, only when the group count equals the expected coordination number), so bent metallocenes (TiCat 14→4, ferrocene 10→2, η²-alkene 5→4) become eligible for selection instead of falling back to the lowest-energy conformer (which sometimes lands a TPY-ish embed). Strictly non-regressive (falls back to lowest-energy unless a conformer both classifies as the target *and* fits tighter); TiCat1/3 now hold `[Ti_TET]` deterministically (8/8 FF, RMSD 0.05–0.23; DEBUG confirms selection actively picks a non-lowest-energy rank), with Ferrocene/TiCp2Me2/Zeise non-regressive.

### Changed
- **`generate_3d_structures`**: Signature gained `uff_pool_size=50, rmsd_threshold=0.5, energy_threshold=2.0`; the conformer pool is now energy-sorted and deduplicated before selection.
- **`TMCOptimizer` (`clean_geometry.py`)**: `clean_geometry()`/`ff_clean()` now return `(success, final_energy)` and stamp `.energy` on the molecule (propagated through `MetalComplex.get_molecule()`), so conformers carry a rankable energy.

## [0.3.1] - 2026-07-05

### Added
- **`MetalloGenAdapter`**: Added support for MACE MLIP optimizer (`mace-omol-0-extra-large-1024`) to refine structures.
- **Force Field configuration**: Surfaced FF convergence knobs via presets, environment variables, and CLI (`--ff-preset`).

### Changed
- **`rmsd_utils.py`**: Replaced generic `999.0` return codes with distinct descending error codes (`998.0`, `997.0`, etc.) for easier debugging.
- **`MetalloGenAdapter`**: Mapped `FF` and `none` (case-insensitive) to `None` for the optimizer flag to default to FF-relaxed geometry.

### Fixed
- **OIN encoder**: Mapped binding atoms to SMILES index via canonical output order, fixing Ir `[cH]` drift.
- **`MetalloGenAdapter`**: Carried encoded sp3-carbon stereo into contract mol, fixing BDPP/BDNN round-trip failures.
- **`MetalloGenAdapter`**: Implemented manual bond-order transfer, fixing CO encoding for FeCO5/FeH2CO4.
- **`MetalloGenAdapter`**: Completed classification audit for MACE geometries, specifically handling SPY pucker and TiCat3/4 TET goldens.

## [0.3.0] - 2026-07-04

### Fixed
- **OIN v3.7: descriptor-free metal token** (`[Pt_SPL]`, was `[Pt@SP1_SPL]`). The `@desc` was an RDKit non-tetrahedral stereo leak via a stale `is_metal` variable in xyz2mol.py; isomer information was and remains fully encoded by slot ordering. Parsers continue to accept legacy `@desc` strings.
- **Eta-ring canonicalization** (`utils/oin_aligner.py`): multi-substituted η-rings now round-trip byte-stably. Fragment order for same-mass haptic ligands is keyed on a heading-independent canonical ring SMILES (was xyz2mol arrival order), and the heading/marker atom of a substituted η-ring not in `SYMMETRIC_LIGANDS` is chosen by lowest `Chem.CanonicalRankAtoms` rank (was 3D geometric alignment to the template slot vector, which varied between a hand-built and a generated structure). Winding-sign computation is untouched; symmetric-η and non-η fragments are unchanged by construction. Dead `base_sort_key` retired.
- **Square-planar Zone-A P enforcement was one-sided** (`generation/molassembler_adapter.py`, `core/chirality.py`): on SPL complexes a metal-bound P stereocenter could only ever generate one enantiomer (the other emitted the wrong 3D and a "could not be enforced" warning), because the metal-present CIP is fixed by placement geometry while the enforcement loop only re-seeded the embed. Fixed by embedding the fragment with a Z=0 dummy metal so `[P@]`/`[P@@]` embed as true 4-coordinate mirror images (symmetric with the encode side). Both enantiomers now generate correctly with no warning.
- **Bidentate incompatible-bite chelates route to the DG fallback** (`generation/molassembler_adapter.py`): DIPAMP-class ligands whose isolated conformation cannot span the chelate bite were placed on the template path and collided with the metal (non-binding H atoms landing ~1.4–1.65 Å from the metal, later misread as hydrides). A non-binding-H proximity guard now routes them to distance geometry, which round-trips them byte-identically.

### Added
- **Winding round-trip preservation (Stereo Phase 1)** (`oin/inline.py`, `generation/oin_parser.py`): the slot-tag parser now captures η-ligand winding markers (`{n>}` CW / `{n<}` CCW) and threads them through to `ParsedOIN.winding_by_slot`, so winding survives XYZ→OIN parsing into the 3D generator instead of being silently dropped.
- **Haptic-face control on 3D generation (Stereo Phase 3)** (`generation/molassembler_adapter.py`, `oin/winding.py`): the winding marker now steers which face of an η-ring the metal binds. A signed-circulation check per ring mirrors the fragment across the ring plane when its embedded winding disagrees with the marker (a proper, CIP-invariant correction).
- **Zone-A P stereocenter encoding (Stereo Phase 4)** (`core/chirality.py`): phosphorus stereocenters bonded directly to the metal are encoded as `[P@]`/`[P@@]` using a lone-pair CIP convention derived from a dummy-metal copy, and verified/enforced on regeneration. (Zone-A **N** encoding remains deferred — RDKit clears trivalent `[N@]` amine tags, so it needs an out-of-band marker.)
- **Direct Parser Fragment Mapping (v0.2.2 Blocker #1) audit completion** (2026-05-10): `_extract_oin_constraints()` audited and verified. Returns 3-tuple `(stripped_smiles, constraints_dict, fragment_to_atom_mapping)` for downstream eta-bond and polydentate-ligand processing. Fragment mapping associates OIN fragment ranks to atom indices in the connected SMILES. Renamed from public `extract_oin_constraints` to private `_extract_oin_constraints` (never a public API; 31 total call sites updated, zero unprefixed references remaining). Verification spike `tools/verify_metal_first.py` confirms metal-first invariant on 6 baseline fixtures; 3 Pd chirality test fixtures documented in `tests/fixtures/_exclusions.yml` (all verified geometrically valid via round-trip RMSD < 1.0 Å). Audit tool `tools/audit_extract_calls.py` confirms rename completeness. 55/55 tests passing (5 new fragment mapping tests verify determinism, cisplatin/polydentate correctness, contiguous atom indices, and metal-at-fragment-zero invariant). Hypergraph node `atom_direct_parser_regex` updated with new output type and status set to `clean`. MiniPRD archived.
- **Direct Parser Molassembler Instantiation audit completion** (2026-05-06): `MiniPRD_DirectParser_MolassemblerInstantiation.md` audited and verified. 20/20 unit tests passing (deterministic Cisplatin/TiCat1 construction, shape assignment, eta bond handling, error cases, all-or-nothing semantics). Implementation in `src/oinsmiles/generation/oin_parser.py` includes `construct_molassembler_mol()` (all-or-nothing transaction wrap), `convert_bond_type()` (RDKit→SCINE mapping), and `extract_oin_constraints()` (OIN v3.6 annotation extraction). SCINE shape mapping covers 10 geometries (SQP, SPL, OCT, TBP, LIN, TPL, TET, TPY, SPY, PBP). New hypergraph node `atom_direct_parser_masm` added to `architecture.yml`. MiniPRD archived.
- **Direct Parser AST Tokenization audit completion** (2026-05-06): `MiniPRD_DirectParser_ASTTokenization.md` audited and verified. 16/16 unit tests passing (deterministic atom/bond extraction, aromatic preservation, implicit H handling, error cases). Implementation in `src/oinsmiles/generation/oin_parser.py::tokenize_unsanitized_smiles()` confirmed to parse unsanitized SMILES with RDKit atom maps, preserve aromatic flags, and defer validation to Molassembler. Hypergraph node `atom_direct_parser_ast` status verified as `clean`. MiniPRD archived.
- **Direct Parser MiniPRD audit completion** (2026-05-06): `MiniPRD_DirectParser_RegexPreprocessor.md` audited and spec-aligned to OIN v3.6 inline format. All 14 unit+integration tests passing. New hypergraph node `atom_direct_parser_regex` added to `architecture.yml`. Updated 4 related MiniPRDs (AST Tokenization, Molassembler Instantiation, Integration, Verification) for consistent constraint dict keys and format examples.
- **MiniPRD audit completion** (all v0.2.0 release specs audited, 2026-05-05): All 5 core feature MiniPRDs now audited and archived — Molassembler Spike, Molassembler Adapter, Chiral Encoding, Chiral Tests, CLI. Updated MiniPRD_MolassemblerAdapter Test 5 to reflect core baseline (5 Pt/Fe/Ir complexes); v0.2.1 eta-ligand regressions documented as known limitations.
- **`GeneratedStructure` dataclass** (`generation/molassembler_adapter.py`, re-exported from `engine.py`): `OIN3DGenerator.generate()` now returns `GeneratedStructure(xyz: str, mol: Optional[Chem.Mol])` instead of a plain string. `mol` carries full RDKit bond connectivity and a 3D conformer with the template-placed positions, enabling callers to write MOL/SDF files with proper bond tables.
- **Bond-preserving MOL/SDF output** in QA test scripts: `verify_roundtrip.py` and `compare_dg_strategies.py` now use `gen_result.mol` for MOL/SDF file output so generated structures include bond connectivity (N–H, M–Cl, M–N dative bonds, etc.). XYZ-only mols are still used for RMSD calculation where matching topology is required.
- **`--include-tmqm` flag** (`verify_xyz_to_oin.py`): tmQM examples are now opt-in. The fast script (`run_verification_fast.sh`) and default roundtrip script exclude the ~103-example tmQM dataset; `run_verification_ALL.sh` includes it.
- **DG strategy comparison script** (`compare_dg_strategies.py`): benchmarks `single`, `ensemble`, and `directed` conformer strategies side-by-side on all curated examples with RMSD, min-distance, and timing metrics. Integrated into all three verification scripts.
- **`Ex{N}_{Name}_` prefixed output files**: verification scripts write named output artifacts (e.g. `Ex1_CisPlatinXYZ-OIN-SMILES_original.xyz`, `…_single.mol`, `…_generated.sdf`) for human QA.
- **Molassembler input diagnostics** in `verify_roundtrip.py` Step 2: logs parsed OIN geometry code, fragment/slot assignments, connected SMILES, permutation index, trans-sym pairs, and expected binding atoms before generation.
- **P/N stereocenter test fixtures** (`tests/integration/`): Added three Pd complex fixtures to verify chirality encoding — PdCl2-R-BINAP (axial-chiral BINAP), PdCl2-RR-BDNN (N-chiral diphosphine), PdCl2-RR-BDPP (P-chiral diphosphine). All pass round-trip verification and extend integration test coverage to 25 examples.

### Fixed
- **`cli.py` `oin2xyz` command**: Fixed `_cmd_oin2xyz` to access `.xyz` attribute from `GeneratedStructure` return value. Prior to fix, the function treated the return value as a plain string, causing `AttributeError` after `OIN3DGenerator.generate()` return type changed in v0.2.0.
- **TiCat1/3/4 3D structure generation** (`generation/molassembler_adapter.py`): `_stitch_multi_eta_fragment` was failing to generate 3D coordinates for ansa-metallocenes with aromatic η5 ligands (Cp, indenyl). Root cause: Phase 4 attempted to kekulize extracted ring SMILES (`[cH]1[cH][cH][cH][cH]1`), which fails for 5-membered all-carbon aromatic rings (5π electrons violates Hückel's 4n+2 rule). Solution: replaced with ETKDG embedding on the **full bridged fragment** (both rings + Si + methyls) with de-aromatization (aromatic bonds→SINGLE, clear aromatic flags). Phase 5 was rewritten to extract ring positions directly from the ETKDG conformer and transform them via centroid/plane alignment. Phase 7 methyl placement corrected: removed spurious `*2.0` scaling and fixed H direction sign (`-cos(tet_angle)` not `cos(tet_angle)`). Result: TiCat1/3/4 now generate 3D with correct atom counts, Si–C bonds (1.87 Å), and tetrahedral methyls. Known trade-off: de-aromatization causes round-trip bonding inference to fail (SINGLE instead of AROMATIC), but geometry quality is good (RMSD ~1.6 Å vs prior ~999 Å failures). See `docs/ETKDG_AROMATIC_FIX.md` for full technical details.

### Changed
- `OIN3DGenerator.generate()` return type changed from `str` to `GeneratedStructure`. Callers that previously used the return value as a string should access `.xyz` for the XYZ block.
- `_stitch_fragment()` and `_stitch_eta_fragment()` now return a 3-tuple `(positions, symbols, mol)` instead of a 2-tuple. `mol` is the RDKit mol with bond topology; for `_stitch_eta_fragment` it is `None` when the analytic geometry fallback is used (e.g. Cp anion ligands in ferrocene).
- `_template_generate()` return type changed from `str | None` to `tuple[str, Chem.Mol | None] | None`. Builds a combined RDKit mol by `CombineMols`-ing the metal atom and each fragment mol, adding dative metal–ligand bonds, and setting a conformer from the final `all_pos` array.

## [0.2.0] - 2026-03-07

### Added
- **SCINE Molassembler backend** (`generation/molassembler_adapter.py`): template-based 3D placement for all ligand types; DG fallback for remaining conformers. Replaces Architector entirely.
- **P/N stereocenter encoding** (`core/chirality.py`): `CIPAssigner` reads the full-TMC 3D conformer (pre-fragmentation) and stores CIP codes on P/N atoms; `ChiralityRecoveryUtility` verifies/corrects chiral tags post-fragmentation; `PseudoAtomStrategy` provides fallback for uncomputable stereocenters.
- **CLI** (`oin-smiles`): two subcommands — `xyz2oin <path>` and `oin2xyz <oin>` — registered as a package entry point.
- **`MolassemblerTimeoutError`** exported from `generation/engine.py`; `OIN3DGenerator` accepts a `timeout` parameter (default 60 s).
- **OIN v3.6 inline format** as canonical output of `XYZToSMILES.convert()` (e.g. `[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}`; the `@SP1` descriptor was a stale-variable bug, not a v3.6 design element — see v3.7 fix above).
- Integration round-trip tests (`verify_roundtrip.py`): unified XYZ → OIN → XYZ → OIN flow with RMSD < 1.0 Å and string-identity checks.
- Unit tests for chirality encoding, Molassembler adapter, regression stability, and axial chirality.

### Changed
- `OIN3DGenerator.generate()` now returns an **XYZ block string** (was an Architector `Molecule` object).
- `XYZToSMILES.convert()` now runs `CIPAssigner.assign_all()` on the full TMC mol before fragmentation.
- `pyproject.toml`: replaced Architector dependency with `scine-molassembler>=2.0.0`; added `oin-smiles` CLI entry point.
- OIN format examples in README updated from V2.4 sidecar to V3.6 inline.

### Removed
- `generation/architector_adapter.py` — Architector integration removed.
- `generation/wrapper.py` — Architector wrapper removed.
- `tests/unit/test_architector.py` — superseded by Molassembler adapter tests.

## [0.1.0] - 2025-xx-xx

Initial release with OIN v2.4 sidecar format, Architector backend, and `XYZToSMILES`/`OIN3DGenerator` APIs.
