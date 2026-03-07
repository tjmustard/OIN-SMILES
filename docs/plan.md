# Implementation Plan

**Current Phase:** Phase 2 (Verification & Fixes)
**Progress:** 70% (Core Logic Implemented)

## Goal
Verify that the `OIN-SMILES` engine correctly canonicalizes 3D structures and reconstructs them via Architector without information loss.

## Phases

### Phase 1: Foundation (Completed)
- [x] Project Structure (uv, src layout)
- [x] Core Data Structures (Vectors, Templates)
- [x] RDKit Integration (`xyz2mol`, Sanitizer)

### Phase 2: Core Logic (Implemented, In Verification)
- [x] `XYZToSMILES` (Canonicalizer)
- [x] `ArchitectorAdapter` (Haptic Expansion)
- [ ] **Verification**: Ensure `Cp` and `Allyl` ligands expand correctly.

### Phase 3: Verification Loop (Next Steps)
- [ ] **Run Round-Trip Tests**: `tests/run_roundtrip.sh`
  - Targets: `Ferrocene.xyz`, `TiCp2Me2.xyz`, `CisPlatin.xyz`.
- [ ] **Fix Regressions**:
  - Address any "SMILES Drift" or "Vector Mismatch" issues.

## Verification Plan
### Automated Tests
- **Command**: `bash tests/run_roundtrip.sh`
- **Success Criteria**:
  - 0% Error Rate on `verify_roundtrip.py`.
  - Generated XYZs must match input Geometry (RMSD < 0.2A).
