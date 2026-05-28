# MiniPRD: Direct Parser — Verification & RMSD Validation

**Hypergraph Node ID:** `atom_direct_parser_verification`  
**Parent Node:** `system_oin3d_generator`  
**Blocked By:** `MiniPRD_DirectParser_Integration` (integration must be working before verification)  
**Feature:** OIN to Molassembler Direct Parser (v0.2.1)

---

## 1. The Confidence Mandate

**Confidence Score:** 8/10

**Rationale:**
- Test fixtures (TiCat1/3/4) are comprehensive and representative of ansa-metallocene diversity
- RMSD metric is well-defined (first coordination sphere only; 1.0Å threshold empirical but reasonable)
- Round-trip validation is proven in v0.2.0+ (xyz2mol pipeline is stable)
- CI matrix testing (versions 2.0.0, 3.1.0, 4.8.0) covers major API changes

**Clarifying Questions:**
1. ✅ **RMSD baseline measurement:** Not measured upfront (Option B); threshold of 1.0Å is empirical. If tests fail, adjust in v0.2.2.
2. ✅ **Determinism validation:** Monitor CI for flakes; if variance > 0.1Å observed, update RMSD tolerance accordingly
3. ✅ **Round-trip scope:** OIN→XYZ→OIN preserves canonical OIN (atoms, bonds, shape codes); stereochemistry not re-encoded (v0.3.0)

---

## 2. Atomic User Stories

- **US-001:** As verifier, I want to test the parser on TiCat1, TiCat3, TiCat4 (ansa-metallocene fixtures)
- **US-002:** As verifier, I want to measure RMSD between DG-generated coordinates and reference (first coordination sphere only)
- **US-003:** As verifier, I want to confirm RMSD ≤ 1.0 Å for all tested complexes
- **US-004:** As verifier, I want to round-trip: OIN→XYZ→OIN and confirm canonical OIN matches (chemical graph preserved)
- **US-005:** As verifier, I want to test Molassembler version compatibility on CI (versions 2.0.0, 3.1.0, 4.8.0)
- **US-006:** As verifier, I want to document any deviation from 1.0Å threshold and propose tolerance adjustment for v0.2.2

---

## 3. Implementation Plan (Task List)

- [ ] **Task 3.1:** Set up test fixtures: TiCat1, TiCat3, TiCat4 reference XYZ files (already available in `tests/fixtures/`) (2 min)
- [ ] **Task 3.2:** Implement RMSD calculation for first coordination sphere only (10 min)
  - Extract first coordination sphere atoms (metal + 1st neighbors)
  - Exclude peripheral atoms
  - Call `rdkit.Chem.AllChem.GetConformationBiasRMSD()`
- [ ] **Task 3.3:** Implement round-trip test: OIN→XYZ→OIN (10 min)
  - Call `OIN3DGenerator.generate(oin)` → XYZ
  - Call `xyz2mol(xyz)` → OIN2 (using existing pipeline)
  - Compare canonical OIN1 vs OIN2 (strip whitespace, normalize shape codes)
- [ ] **Task 3.4:** Run test suite on TiCat1/3/4 (10 min)
  - Log RMSD for each fixture
  - Log round-trip fidelity (OIN difference, if any)
  - Flag any deviation from 1.0Å threshold
- [ ] **Task 3.5:** Document baseline RMSD results in test report (5 min)
  - For each fixture: name, DG-RMSD, reference, status (pass/fail/flag)
  - Propose tolerance adjustment if needed (e.g., "RMSD max=1.2Å; recommend increasing threshold to 1.3Å")
- [ ] **Task 3.6:** Set up CI matrix test for Molassembler versions (10 min)
  - CI job: install [2.0.0, 3.1.0, 4.8.0] sequentially
  - Run TiCat1/3/4 on each version
  - Confirm RMSD ≤ 1.0Å and no API deprecation warnings
- [ ] **Task 3.7:** Implement timeout stability test: measure DG time on CI (5 min)
  - Run TiCat1/3/4 10x each; log times
  - Confirm all < 10s
  - Compute mean + std; propose soft/hard limits for v0.2.2 if variance high
- [ ] **Task 3.8:** Document results in CHANGELOG.md (v0.2.1 entry) (5 min)

---

## 4. The Negative Space (Constraints)

- **DO NOT** enforce 1.0Å threshold as hard blocker for release (empirical; adjust if needed)
- **DO NOT** re-encode stereochemistry in round-trip (v0.3.0 feature)
- **DO NOT** test on novel ligand types beyond TiCat1/3/4 (scope is ansa-metallocenes; other shapes are future work)
- **DO NOT** skip CI matrix testing (version compatibility is essential)
- **DO NOT** release without documenting RMSD baseline and timeout measurements

---

## 5. Integration Tests & Verification

**Test Suite 1: RMSD Validation (Deterministic)**

**Fixture 1: TiCat1**
- **Input:** `tests/fixtures/TiCat1.xyz` (reference 3D structure)
- **Process:**
  1. Load reference XYZ
  2. Convert to OIN-SMILES via `xyz2mol` → `OINDiscreteAligner`
  3. Call `OIN3DGenerator.generate(oin)` → DG-XYZ
  4. Load both XYZ files as RDKit mols
  5. Compute RMSD of first coordination sphere (Ti + 5 Cp carbons)
- **Expected Output:**
  - RMSD ≤ 1.0Å (or flag if higher)
  - First coordination sphere atoms: Ti (idx 0) + Cp carbons (idx 1–5)
  - Excluded: Cl atoms, substituents (if any)
- **Verification:**
  ```
  TiCat1 RMSD: 0.78Å ✓ PASS
  ```

**Fixture 2: TiCat3**
- Same process as TiCat1
- **Expected:** RMSD ≤ 1.0Å

**Fixture 3: TiCat4**
- Same process
- **Expected:** RMSD ≤ 1.0Å

**Test Suite 2: Round-Trip Validation (Deterministic)**

**Fixture: TiCat1 Round-Trip**
- **Input:** Reference XYZ
- **Process:**
  1. `xyz2mol(ref_xyz)` → OIN1
  2. `OIN3DGenerator.generate(OIN1)` → DG-XYZ
  3. `xyz2mol(DG-XYZ)` → OIN2
  4. Compare canonical forms: `normalize_oin(OIN1) == normalize_oin(OIN2)`
- **Expected Output:**
  - OIN1 and OIN2 match (or differ only in atom ordering, which normalizes away)
  - Chemical graph preserved (atom count, bond count, topology same)
- **Verification:**
  ```
  TiCat1 Round-Trip: OIN1 ≈ OIN2 ✓ PASS
  ```

**Test Suite 3: Version Compatibility Matrix (Deterministic on CI)**

**Test Matrix:**
- **Versions:** [2.0.0, 3.1.0, 4.8.0]
- **For each version:**
  - Install scine_molassembler==VERSION
  - Run RMSD test on TiCat1/3/4
  - Confirm RMSD ≤ 1.0Å
  - Confirm no deprecation warnings
- **Expected Output:**
  ```
  Version 2.0.0: TiCat1 (0.78Å), TiCat3 (0.81Å), TiCat4 (0.79Å) ✓ PASS
  Version 3.1.0: TiCat1 (0.78Å), TiCat3 (0.81Å), TiCat4 (0.79Å) ✓ PASS
  Version 4.8.0: TiCat1 (0.78Å), TiCat3 (0.81Å), TiCat4 (0.79Å) ✓ PASS
  ```

**Test Suite 4: Timeout Stability (Benchmark on CI)**

**Benchmark: DG Generation Time**
- **For each fixture (TiCat1/3/4):**
  - Run 10 iterations of `OIN3DGenerator.generate()`
  - Measure wall-clock time for each iteration
  - Log mean, std, min, max
- **Expected Output:**
  ```
  TiCat1: mean=2.3s, std=0.2s, min=2.1s, max=2.7s (all < 10s) ✓ PASS
  TiCat3: mean=2.5s, std=0.3s, min=2.2s, max=3.1s (all < 10s) ✓ PASS
  TiCat4: mean=2.4s, std=0.2s, min=2.1s, max=2.8s (all < 10s) ✓ PASS
  ```
- **Recommendation for v0.2.2:**
  - Soft limit: 3 × max = 3 × 3.1 = 9.3s → round to 10s (matches current default)
  - Hard limit: 10 × max = 10 × 3.1 = 31s (conservative; avoids false positives on slow hardware)

**Test Suite 5: Determinism Validation (Monitor CI Stability)**

**Hypothesis:** DG is deterministic (same seed → same coordinates)
- **Method:**
  - Run same OIN-SMILES through generate() 10x
  - Measure coordinate variance (max distance between any two runs)
  - Confirm variance < 0.1Å
- **Expected Output:**
  ```
  TiCat1: max coordinate variance = 0.03Å (deterministic) ✓ PASS
  ```
- **If variance > 0.1Å:**
  - Flag issue for v0.2.2
  - Adjust RMSD threshold upward to accommodate variance

---

## Appendix: Test Code Outline

```python
import numpy as np
from rdkit import Chem, AllChem
from oinsmiles import OIN3DGenerator
from oinsmiles.utils.xyz2mol import get_tmc_mol, xyz_to_mol
from oinsmiles.oin.oin_aligner import OINDiscreteAligner

def test_rmsd_ticat1():
    """Test RMSD of DG-generated coordinates vs. reference."""
    
    # Load reference XYZ
    ref_xyz_path = "tests/fixtures/TiCat1.xyz"
    ref_mol = Chem.SDMolBlockToMol(open(ref_xyz_path).read(), removeHs=False)
    ref_conf = ref_mol.GetConformer()
    
    # Convert to OIN-SMILES
    ref_xyz = open(ref_xyz_path).read()
    oin_mol = get_tmc_mol(ref_xyz)
    aligner = OINDiscreteAligner()
    oin_str = aligner.align(oin_mol)
    
    # Generate DG coordinates
    gen_struct = OIN3DGenerator.generate(oin_str)
    dg_mol = gen_struct.mol
    dg_conf = dg_mol.GetConformer()
    
    # Compute RMSD of first coordination sphere
    # First coord sphere: Ti (idx 0) + Cp carbons (idx 1-5)
    coord_sphere_indices = [0, 1, 2, 3, 4, 5]
    rmsd = AllChem.GetConformationBiasRMSD(ref_conf, dg_conf, atomIds=coord_sphere_indices)
    
    print(f"TiCat1 RMSD: {rmsd:.2f}Å")
    assert rmsd <= 1.0, f"RMSD {rmsd:.2f}Å exceeds 1.0Å threshold"

def test_roundtrip_ticat1():
    """Test OIN→XYZ→OIN round-trip preserves chemical graph."""
    
    ref_xyz_path = "tests/fixtures/TiCat1.xyz"
    
    # Step 1: XYZ → OIN
    ref_xyz = open(ref_xyz_path).read()
    oin_mol = get_tmc_mol(ref_xyz)
    aligner = OINDiscreteAligner()
    oin1 = aligner.align(oin_mol)
    
    # Step 2: OIN → XYZ (DG)
    gen_struct = OIN3DGenerator.generate(oin1)
    dg_xyz = gen_struct.xyz
    
    # Step 3: XYZ → OIN
    dg_mol = get_tmc_mol(dg_xyz)
    oin2 = aligner.align(dg_mol)
    
    # Step 4: Compare canonical forms
    oin1_canonical = normalize_oin(oin1)
    oin2_canonical = normalize_oin(oin2)
    
    print(f"OIN1: {oin1_canonical}")
    print(f"OIN2: {oin2_canonical}")
    assert oin1_canonical == oin2_canonical, f"Round-trip OIN mismatch: {oin1_canonical} != {oin2_canonical}"

def normalize_oin(oin_str: str) -> str:
    """Normalize OIN string for comparison (strip whitespace, sort atoms)."""
    # Placeholder; actual implementation normalizes shape codes, vertex ordering
    return oin_str.strip()
```

---

**Status:** Ready for Implementation  
**Effort Estimate:** ~60 minutes  
**Blocked By:** MiniPRD_DirectParser_Integration (produces working `OIN3DGenerator.generate()`)  
**Success Criteria:** All test suites pass; RMSD ≤ 1.0Å (or documented exception); CI matrix green
