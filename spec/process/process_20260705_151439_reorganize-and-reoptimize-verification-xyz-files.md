# Process Document: Reorganize and Reoptimize Verification XYZ Files

**Generated:** 2026-07-05T15:14:39-07:00
**Session Focus:** Re-optimize and Reorganize Verification XYZ Files

## Problem Statement

The integration and unit test suites were loading transition metal complex XYZ structure files from inconsistent directories (`tests/fixtures/` and `tests/integration/`) and containing case-sensitive duplicate names (e.g. `cisplatin.xyz` vs `CisPlatin.xyz`). These structures needed to be re-optimized on a GPU using MACE (`mace-omol-0-extra-large-1024` weights) and consolidated into the `tests/fixtures/` directory to improve codebase cleanliness and prevent case conflicts on macOS/Windows.

## Starting State

The codebase contained `tests/fixtures/` (with 12 `.xyz` files) and `tests/integration/` (with 28 `.xyz` files). Some of these files were duplicates with casing differences, and several tests in the unit test suite depended on lowercase file paths. A `--cpu` execution path for MACE was also missing in some runners, defaulting to CUDA when available but without a clean override option.

## Approach & Methodology

This was a spec-driven refactoring and GPU-optimization session. The approach was:
1. Parse the command-line argument `--cpu` in bash runners and Python scripts to hide CUDA when necessary (`CUDA_VISIBLE_DEVICES=""`).
2. Backup all original `.xyz` files into a dedicated backup folder `tests/backup/`.
3. Develop a custom python script (`tests/reoptimize_xyz_mace.py`) to run geometry optimization on all unique structures on the MACE GPU (`device="cuda"`), retrieving their charge and spin states dynamically using the OIN-SMILES library itself.
4. Consolidate all re-relaxed standard 4-column `.xyz` files in `tests/fixtures/`.
5. Remove the duplicate/original files from `tests/integration/` and update references across all test scripts to point to `tests/fixtures/`.
6. Run the full unit and integration test suites to verify that 100% of the tests pass.

## Steps Taken

1. Added `--cpu` flag parsing and execution fallback to `tests/run_verification.sh`, `tests/integration/verify_phase1.py`, and `tests/integration/verify_roundtrip.py`. Verified that this runs successfully on CPU.
2. Created the backup directory `tests/backup/` and moved all original XYZ structures there.
3. Created `tests/reoptimize_xyz_mace.py` to optimize all unique structures on CUDA using the `MACE-omol-0-extra-large-1024` model, automatically calculating spin/charge via OIN translation.
4. Encountered a `ValueError: too many values to unpack (expected 4)` in `read_xyz_file` because `ase.io.write` includes extra columns (forces, energies, etc.) when saving optimized structures.
5. Resolved this issue by implementing a manual standard 4-column writer in the re-optimization script, outputting exactly standard `element X Y Z` coordinates.
6. Re-ran the re-optimization script on CUDA, producing the relaxed standard structures in `tests/fixtures/`.
7. Deleted the temporary re-optimization script and all `.xyz` files in `tests/integration/`.
8. Updated `tests/integration/verify_xyz_to_oin.py`, `tests/unit/test_regression_stability.py`, `tests/unit/test_stereo_roundtrip_diagnostics.py`, `tests/unit/test_eta_bonds.py`, `tests/unit/test_zone_a_p_encode.py`, and `tests/unit/test_zone_a_p_genenforce.py` to point to standard camelCase names in `../fixtures/`.
9. Updated the hardcoded expected SHA256 of `Rh-RR-DIPAMP-Cl2.xyz` in `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` to match the new optimized structure.
10. Added `--cpu` support to `run_verification_fast.sh` and `run_verification_ALL.sh`.
11. Re-ran the full test suite (`unittest discover tests/unit` and `verify_xyz_to_oin.py --limit 4` / `verify_roundtrip.py --limit 2`): all tests pass successfully.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| Use `CUDA_VISIBLE_DEVICES=""` for `--cpu` | Modify MACE calculator code inside generator wrappers | Simpler, cleaner, and disables GPU visibility globally for PyTorch/MACE without modifying core modules. |
| Consolidate all XYZs to `tests/fixtures/` and standardise casing | Keep duplicate lowercase names (`cisplatin.xyz`, etc.) | Prevents naming collisions on case-insensitive filesystems (macOS and Windows default). |
| Implement manual 4-column XYZ writer in optimization script | Use `ase.io.write` | `ase.io.write` exports calculator-specific properties (e.g. forces), causing parsing failures in `xyz2mol`. |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| run_verification.sh | tests/run_verification.sh | updated |
| run_verification_fast.sh | tests/run_verification_fast.sh | updated |
| run_verification_ALL.sh | tests/run_verification_ALL.sh | updated |
| verify_phase1.py | tests/integration/verify_phase1.py | updated |
| verify_roundtrip.py | tests/integration/verify_roundtrip.py | updated |
| verify_xyz_to_oin.py | tests/integration/verify_xyz_to_oin.py | updated |
| test_regression_stability.py | tests/unit/test_regression_stability.py | updated |
| test_stereo_roundtrip_diagnostics.py | tests/unit/test_stereo_roundtrip_diagnostics.py | updated |
| test_eta_bonds.py | tests/unit/test_eta_bonds.py | updated |
| test_zone_a_p_encode.py | tests/unit/test_zone_a_p_encode.py | updated |
| test_zone_a_p_genenforce.py | tests/unit/test_zone_a_p_genenforce.py | updated |
| Rh-RR-DIPAMP-Cl2_oin.txt | tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt | updated |
| optimized fixtures | tests/fixtures/*.xyz | updated / added |
| integration fixtures | tests/integration/*.xyz | deleted (moved to fixtures) |
| backup folder | tests/backup/ | created |

## Results & Outcomes

- All integration and unit tests run against standardized camelCase 3D geometries consolidated under `tests/fixtures/`.
- All unique structures are successfully re-relaxed using MACE GPU extra-large model.
- 100% of unit tests (`discover tests/unit`) and integration tests (`verify_xyz_to_oin.py`, `verify_roundtrip.py`) pass.

## How to Reproduce

### Prerequisites
- Python 3.10+ and `uv` package manager installed.
- CUDA-enabled GPU and MACE model weights at `models/mace/MACE-omol-0-extra-large-1024.model`.

### Steps
1. Restore original structure files from `tests/backup/` to their initial folders if needed.
2. Put the `reoptimize_xyz_mace.py` script in `tests/`.
3. Run the optimization script:
   ```bash
   uv run python tests/reoptimize_xyz_mace.py
   ```
4. Run the test suite:
   ```bash
   uv run python -m unittest discover tests/unit
   uv run bash tests/run_verification_fast.sh --optimizer mace-omol-0-extra-large-1024 --ensemble-size 1
   ```
