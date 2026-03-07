# Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

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
