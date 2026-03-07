# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-03-07

### Added
- **CLI entry point**: `oin-smiles xyz2oin <path>` converts XYZ to OIN-SMILES; `oin-smiles oin2xyz <oin>` generates XYZ from OIN. Non-zero exit codes on errors.
- **`CIPAssigner`**: Pre-fragmentation CIP assignment on full sanitized mol. Stores `_CIPCode` on P/N atoms before fragmentation destroys metal context.
- **`ChiralityRecoveryUtility`**: Re-applies `@`/`@@` chirality tags to ligand SMILES after `OINSanitizer` using stored `_CIPCode` properties.
- **`PseudoAtomStrategy`**: Lone-pair pseudo-atom fallback for 3-coordinate P/N stereocenters where RDKit CIP assignment fails. Pseudo-atoms stripped from final OIN output.
- **Axial chirality (`@b` tag)**: Manual biaryl detection via dihedral angle assigns `STEREOATROP_CW/CCW`; serialized as `@b u-v:STEREO` sidecar in OIN string.
- **`MolassemblerAdapter`**: SCINE Molassembler backend for 3D generation. Replaces Architector. Runs in `ProcessPoolExecutor` for GIL-safe 60s timeout.
- **`MolassemblerTimeoutError`**: Raised on DG timeout; triggers exit code 2 from CLI.
- **Chiral test suite**: `test_regression_stability.py`, `test_binap_stability.py`, `test_chiral_p.py`, `test_chiral_n.py`, `test_axial_chiral.py`.
- **`tests/test_helpers.py`**: `extract_ligand_smiles()` test utility for CIP oracle validation.

### Changed
- **`OIN3DGenerator`**: Rewired to use `MolassemblerAdapter` instead of Architector.
- **`OINInlineHandler.parse_inline_string()`**: Regex-only slot marker stripping — no RDKit round-trip to preserve `@`/`@@`.
- **`XYZToSMILES.convert()`**: Now calls `CIPAssigner.assign_all()` on full mol before fragmentation.
- **`Canonicalizer.canonicalize()`**: Collects `OIN_BondStereo` from bonds and emits `@b` sidecar tag.
- **`OINParser`**: Added `parse_bond_stereo()` to parse the new `@b` tag.
- **`OINInlineHandler`**: Preserves `@b` tag in `generate_inline_string` output.
- **`pyproject.toml`**: Replaced `architector` with `scine-molassembler>=3.0.0`. Added `[project.scripts]` entrypoint.
- **`tests/`**: Removed obsolete `test_adapter.py`, `test_architector.py`, `test_roundtrip.py` (referenced deleted Architector code).

### Removed
- **`ArchitectorAdapter`** (`src/oinsmiles/generation/architector_adapter.py`) — deleted.
- **`ArchitectorWrapper`** — deleted.
- All `architector` and `xtb` runtime dependencies removed from generation path.

## [0.1.0] - 2026-02-16

### Added
- Initial release: OIN-SMILES V3.6 Inline Syntax.
- XYZ → OIN-SMILES pipeline with PAI alignment, fragmentation, and haptic vector encoding.
- OIN-SMILES → 3D pipeline via Architector backend.
- Geometry support: `SPL`, `OCT`, `TET`, `TBP`, `LIN`, `TPY`, `SPY`, `TPL`, `PBP`.
- Haptic directionality tags (`<`, `>`).
- Zone A sanitization for coordinating atoms.
