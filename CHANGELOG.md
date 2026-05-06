# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- **`GeneratedStructure` dataclass** (`generation/molassembler_adapter.py`, re-exported from `engine.py`): `OIN3DGenerator.generate()` now returns `GeneratedStructure(xyz: str, mol: Optional[Chem.Mol])` instead of a plain string. `mol` carries full RDKit bond connectivity and a 3D conformer with the template-placed positions, enabling callers to write MOL/SDF files with proper bond tables.
- **Bond-preserving MOL/SDF output** in QA test scripts: `verify_roundtrip.py` and `compare_dg_strategies.py` now use `gen_result.mol` for MOL/SDF file output so generated structures include bond connectivity (N–H, M–Cl, M–N dative bonds, etc.). XYZ-only mols are still used for RMSD calculation where matching topology is required.
- **`--include-tmqm` flag** (`verify_xyz_to_oin.py`): tmQM examples are now opt-in. The fast script (`run_verification_fast.sh`) and default roundtrip script exclude the ~103-example tmQM dataset; `run_verification_ALL.sh` includes it.
- **DG strategy comparison script** (`compare_dg_strategies.py`): benchmarks `single`, `ensemble`, and `directed` conformer strategies side-by-side on all curated examples with RMSD, min-distance, and timing metrics. Integrated into all three verification scripts.
- **`Ex{N}_{Name}_` prefixed output files**: verification scripts write named output artifacts (e.g. `Ex1_CisPlatinXYZ-OIN-SMILES_original.xyz`, `…_single.mol`, `…_generated.sdf`) for human QA.
- **Molassembler input diagnostics** in `verify_roundtrip.py` Step 2: logs parsed OIN geometry code, fragment/slot assignments, connected SMILES, permutation index, trans-sym pairs, and expected binding atoms before generation.
- **P/N stereocenter test fixtures** (`tests/integration/`): Added three Pd complex fixtures to verify chirality encoding — PdCl2-R-BINAP (axial-chiral BINAP), PdCl2-RR-BDNN (N-chiral diphosphine), PdCl2-RR-BDPP (P-chiral diphosphine). All pass round-trip verification and extend integration test coverage to 25 examples.

### Fixed
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
- **OIN v3.6 inline format** as canonical output of `XYZToSMILES.convert()` (e.g. `[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}`).
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
