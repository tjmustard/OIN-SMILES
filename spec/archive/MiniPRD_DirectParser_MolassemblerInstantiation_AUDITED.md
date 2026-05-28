# MiniPRD: Direct Parser — Molassembler Instantiation

**Hypergraph Node ID:** `atom_direct_parser_masm`  
**Parent Node:** `system_oin3d_generator`  
**Blocked By:** `MiniPRD_DirectParser_ASTTokenization` (AST output is input to Molassembler stage)  
**Feature:** OIN to Molassembler Direct Parser (v0.2.1)

---

## 1. The Confidence Mandate

**Confidence Score:** 8/10

**Rationale:**
- Molassembler API is proven in production (v0.2.0+, MolassemblerAdapter)
- All-or-nothing transaction semantics are achievable with try/finally
- Error messages from Molassembler are opaque; translation required but straightforward

**Clarifying Questions:**
1. ✅ **Timeout handling during DG generation:** Default 10s (resolve in Integration MiniPRD); hotfix v0.2.2
2. ✅ **Determinism of DG engine:** Assume deterministic; validate in Verification MiniPRD via CI stability

---

## 2. Atomic User Stories

- **US-001:** As Molassembler instantiator, I want to create an empty `scine_molassembler.Molecule()` object
- **US-002:** As instantiator, I want to add atoms sequentially via `mol.push_back(atom)` and preserve index ordering
- **US-003:** As instantiator, I want to add standard bonds (single, double, aromatic) via `mol.add_bond(i, j, bond_type)`
- **US-004:** As instantiator, I want to assign polyhedral shape to metal center via `mol.set_shape(metal_idx, SCINE_SHAPE)`
- **US-005:** As instantiator, I want to add eta bonds from metal to ligands via `mol.add_bond(metal_idx, lig_idx, BondType.eta)`
- **US-006:** As error handler, I want to detect construction failures and raise user-friendly `ValueError` with context
- **US-007:** As transaction manager, I want to ensure all-or-nothing behavior: if any step fails, the entire molecule is discarded

---

## 3. Implementation Plan (Task List)

- [ ] **Task 3.1:** Define SCINE shape library mapping (`SQP` → `SCINE_SHAPES.SquarePlanar`, etc.) (5 min)
- [ ] **Task 3.2:** Create empty Molassembler molecule: `mol = masm.Molecule()` (2 min)
- [ ] **Task 3.3:** Implement atom addition loop: for each atom from AST, call `mol.push_back(atom)` (5 min)
  - Verify index preservation: `mol.NumAtoms()` matches input count
  - Handle atom properties (atomic number, formal charge)
- [ ] **Task 3.4:** Implement standard bond addition: for each bond from AST, call `mol.add_bond(i, j, bond_type)` (5 min)
  - Map RDKit `BondType` to Molassembler `BondType`
  - Validate endpoints exist before adding
- [ ] **Task 3.5:** Implement shape assignment: look up metal center in constraints, call `mol.set_shape(metal_idx, SCINE_SHAPE)` (5 min)
  - Error handling: if shape not in library, catch and re-raise with list of valid shapes
- [ ] **Task 3.6:** Implement eta bond addition: for each eta index in constraints, call `mol.add_bond(metal_idx, lig_idx, BondType.eta)` (5 min)
- [ ] **Task 3.7:** Implement all-or-nothing error handling: try/finally or copy-construct on failure (10 min)
- [ ] **Task 3.8:** Unit test: construct molecules for Cisplatin, TiCat1/3/4; verify atom counts, bond types, shape assignments (15 min)

---

## 4. The Negative Space (Constraints)

- **DO NOT** sanitize the RDKit atoms before adding to Molassembler
- **DO NOT** validate bond orders or valence (Molassembler will reject invalid configuration)
- **DO NOT** assume metal center is at atom index 0 (look up in constraints dict)
- **DO NOT** return half-constructed molecule on error; use all-or-nothing semantics
- **DO NOT** attempt fallback strategies if shape assignment fails; fail fast with error message

---

## 5. Integration Tests & Verification

**Test 1 (Deterministic — Cisplatin):**
- **Input (from ASTTokenization):**
  ```python
  atoms = [Pt, Cl, Cl, N, N]  # from RDKit
  bonds = [(0,1), (0,2), (0,3), (0,4)]  # Pt-Cl, Pt-Cl, Pt-N, Pt-N
  constraints = {0: {'shape': 'SQP', 'vertex_indices': [1, 2, 3, 4]}}
  ```
- **Expected Output:**
  ```python
  mol = scine_molassembler.Molecule()
  mol.NumAtoms() == 5
  mol.NumBonds() == 4
  # Shape: metal (idx 0) has SquarePlanar geometry
  ```
- **Verification:**
  - 5 atoms added without error
  - 4 bonds added (all single)
  - Shape assigned to metal center (idx 0)

**Test 2 (Deterministic — TiCat1 with Eta Bonds):**
- **Input:**
  ```python
  atoms = [Ti, Cp (5 atoms), Cl, Cl]  # 8 atoms total
  bonds = [(0,1), (0,2), (0,3), (0,4), (0,5), (1,2), (2,3), (3,4), (4,5), (5,1)]  # Ti-Cp bonds + Cp ring
  constraints = {0: {'shape': 'OC', 'vertex_indices': [1, 2, 3, 4, 5]}}
  ```
- **Expected Output:**
  - 8 atoms added
  - 10 bonds added (Pt-Cp single bonds + Cp ring bonds)
  - 5 eta bonds added (metal to each Cp carbon)
  - Shape: metal (idx 0) has Octahedral geometry
- **Verification:**
  ```python
  mol.NumAtoms() == 8
  mol.NumBonds() == 15  # 10 input + 5 eta
  ```

**Test 3 (Edge Case — Shape Not in Library):**
- **Input:**
  ```python
  constraints = {0: {'shape': 'INVALID_SHAPE', 'vertex_indices': []}}
  ```
- **Expected Behavior:**
  - `mol.set_shape(0, SCINE_SHAPES.INVALID_SHAPE)` raises exception
  - Exception caught, re-raised as:
    ```python
    ValueError("Unknown shape INVALID_SHAPE; valid shapes: [SQP, OC, TBP, ...]")
    ```

**Test 4 (Error Case — Bond Endpoint Out of Bounds):**
- **Input:**
  ```python
  atoms = [Pt, Cl]  # 2 atoms
  bonds = [(0, 1), (0, 99)]  # Second bond references atom 99 (doesn't exist)
  ```
- **Expected Behavior:**
  - During bond addition, Molassembler validation fails
  - Exception caught, re-raised as:
    ```python
    ValueError("Failed to add bond (0, 99): atom index 99 out of bounds (mol has 2 atoms)")
    ```

**Test 5 (All-or-Nothing Transaction):**
- **Input (simulated failure at step 5/7):**
  - Add atoms 0–4 successfully
  - Add bonds 0–3 successfully
  - Try to add bond (0, 99) → **FAILS**
- **Expected Behavior:**
  - Exception caught in try/except
  - Molecule discarded (not returned)
  - User gets `ValueError` with context about which step failed
  - **Verification:** No partial molecule returned; callers get exception or valid mol, never half-constructed state

---

## Appendix: Code Outline

```python
import scine_molassembler as masm
from typing import List, Tuple, Dict

SCINE_SHAPE_MAP = {
    'SQP': masm.shapes.Shapes.SquarePlanar,
    'OC': masm.shapes.Shapes.Octahedral,
    'TBP': masm.shapes.Shapes.TrigonalBipyramidal,
    # ... etc
}

def construct_molassembler_mol(
    atoms: List,
    bonds: List[Tuple],
    constraints: Dict,
    mol_rdkit
) -> masm.Molecule:
    """
    Construct Molassembler molecule with atoms, bonds, shape, and eta bonds.
    
    All-or-nothing: if any step fails, raise ValueError and do not return partial mol.
    
    Args:
        atoms: list of RDKit Atom objects
        bonds: list of (i, j, bond_type) tuples
        constraints: {atom_idx: {'shape': str, 'eta_indices': list[int]}}
        mol_rdkit: RDKit molecule (for properties, returned separately)
    
    Returns:
        masm.Molecule with all atoms, bonds, shape, eta bonds added
    
    Raises:
        ValueError: with context about which step failed
    """
    
    try:
        # Step 1: Create empty molecule
        mol = masm.Molecule()
        
        # Step 2: Add atoms
        for atom in atoms:
            mol.push_back(masm.AtomCollection.Atom(
                atomic_number=atom.GetAtomicNum(),
                formal_charge=atom.GetFormalCharge()
            ))
        
        if mol.NumAtoms() != len(atoms):
            raise ValueError(f"Atom count mismatch: added {len(atoms)}, mol has {mol.NumAtoms()}")
        
        # Step 3: Add standard bonds
        for i, j, bond_type in bonds:
            mol.add_bond(i, j, convert_bond_type(bond_type))
        
        # Step 4: Assign shape to metal center
        if 0 in constraints and 'shape' in constraints[0]:
            shape_code = constraints[0]['shape']
            if shape_code not in SCINE_SHAPE_MAP:
                raise ValueError(f"Unknown shape {shape_code}; valid: {list(SCINE_SHAPE_MAP.keys())}")
            scine_shape = SCINE_SHAPE_MAP[shape_code]
            mol.set_shape(0, scine_shape)  # Assume metal is at index 0
        
        # Step 5: Add eta bonds
        if 0 in constraints and 'eta_indices' in constraints[0]:
            for lig_idx in constraints[0]['eta_indices']:
                mol.add_bond(0, lig_idx, masm.BondType.eta)
        
        return mol
    
    except Exception as e:
        # All-or-nothing: raise with context
        raise ValueError(f"Failed to construct Molassembler molecule: {e}") from e

def convert_bond_type(rdkit_bond_type):
    """Map RDKit BondType to Molassembler BondType."""
    # Implementation: RDKit.Single → masm.BondType.Single, etc.
    pass
```

---

**Status:** Ready for Implementation  
**Effort Estimate:** ~55 minutes  
**Blocked By:** MiniPRD_DirectParser_ASTTokenization (produces atoms, bonds)
