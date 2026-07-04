# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

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
