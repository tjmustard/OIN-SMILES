# MiniPRD: Direct Parser — AST Tokenization

**Hypergraph Node ID:** `atom_direct_parser_ast`  
**Parent Node:** `system_oin3d_generator`  
**Blocked By:** `MiniPRD_DirectParser_RegexPreprocessor` (regex output is input to AST stage)  
**Feature:** OIN to Molassembler Direct Parser (v0.2.1)

---

## 1. The Confidence Mandate

**Confidence Score:** 9/10

**Rationale:**
- RDKit's unsanitized SMILES parsing is well-documented and stable
- Implicit hydrogen handling is predictable (RDKit preserves implicit H even without sanitization)
- No aromaticity enforcement at this stage means no surprises
- Pattern: `Chem.MolFromSmiles(smiles, sanitize=False)` is standard in RDKit workflows

**Clarifying Questions:**
None. RDKit API is well-defined; only requirement is to preserve unsanitized state.

---

## 2. Atomic User Stories

- **US-001:** As AST tokenizer, I want to parse unsanitized SMILES strings (with RDKit atom maps) into RDKit `Mol` objects
- **US-002:** As tokenizer, I want to extract atoms and bonds from the RDKit `Mol` in stable order (preserving indices)
- **US-003:** As tokenizer, I want to document and verify implicit hydrogen handling (preserved even without sanitization)
- **US-004:** As tokenizer, I want to ensure aromatic flags are NOT enforced or modified (aromaticity check deferred to Molassembler)
- **US-005:** As error handler, I want to fail gracefully if RDKit rejects the unsanitized SMILES (malformed input)

---

## 3. Implementation Plan (Task List)

- [ ] **Task 3.1:** Call `Chem.MolFromSmiles(stripped_smiles, sanitize=False)` on regex output (2 min)
- [ ] **Task 3.2:** Extract atom list: `mol.GetAtoms()` → iterate, collect atomic number, formal charge, implicit H (5 min)
- [ ] **Task 3.3:** Extract bond list: `mol.GetBonds()` → iterate, collect (i, j, bond_type) tuples (5 min)
- [ ] **Task 3.4:** Verify indices are stable: atom indices match `atom.GetIdx()`, bond endpoints match mol atom count (5 min)
- [ ] **Task 3.5:** Document implicit hydrogen semantics: verify `atom.GetTotalNumHs(implicitH=True)` reflects OIN input (5 min)
- [ ] **Task 3.6:** Unit test: parse 10 OIN-SMILES variants (with and without aromatic flags) and verify bonds/atoms match expected (10 min)
  - Verify aromatic rings are parsed but not kekulized
  - Verify implicit H counts are preserved
  - Verify atom maps from regex stage are intact

---

## 4. The Negative Space (Constraints)

- **DO NOT** sanitize the RDKit molecule (preserve unsanitized state for Molassembler)
- **DO NOT** kekulize aromatic rings (SCINE will handle aromatics in its own format)
- **DO NOT** assign stereochemistry or chirality at this stage (deferred to Molassembler shape assignment)
- **DO NOT** validate bond orders or valence (RDKit's unsanitized parser is lenient; Molassembler will validate)

---

## 5. Integration Tests & Verification

**Test 1 (Deterministic — Simple Complex):**
- **Input (from Regex):** `[Pt:1].[Cl:2].[Cl:3].[N:4].[N:5]`
- **Expected Output:**
  ```python
  atoms = [
    Atom(atomic_num=78, idx=0, formal_charge=0),  # Pt
    Atom(atomic_num=17, idx=1, formal_charge=0),  # Cl
    Atom(atomic_num=17, idx=2, formal_charge=0),  # Cl
    Atom(atomic_num=7, idx=3, formal_charge=0),   # N
    Atom(atomic_num=7, idx=4, formal_charge=0),   # N
  ]
  bonds = [
    Bond(0, 1, BondType.Single),  # If present in SMILES
    # ... etc
  ]
  ```
  - **Verification:** 5 atoms extracted, indices 0–4, implicit H counts preserved

**Test 2 (Deterministic — Aromatic Ring, Unsanitized):**
- **Input (from Regex):** `[Ti:1].[c:2]1[c:3][c:4][c:5][c:6]1`  (aromatic Cp ring, unsanitized)
- **Expected Output:**
  - Atoms extracted: 6 total (Ti + 5 carbons)
  - Bonds extracted: Ti-C bonds + C-C aromatic bonds (marked `BondType.Aromatic`)
  - **Verification:** Aromatic flags present but not kekulized; ring structure intact

**Test 3 (Edge Case — Implicit Hydrogens):**
- **Input:** `[C:1].[C:2]`  (two carbons, no explicit H)
- **Expected Output:**
  ```python
  atoms[0].GetTotalNumHs() == 4  # Implicit H preserved
  atoms[1].GetTotalNumHs() == 4
  ```
  - **Verification:** Implicit hydrogens are preserved even in unsanitized state

**Test 4 (Error Case — Malformed SMILES):**
- **Input:** `[Pt:1].[Cl:99]` (invalid atom map or syntax)
- **Expected Output:** `RDKitException` or `None` from `Chem.MolFromSmiles()`
  - **Handling:** Catch exception, re-raise as `ValueError` with message (deferred to Integration MiniPRD)

---

## Appendix: Code Outline

```python
from rdkit import Chem
from typing import Tuple, List

def tokenize_unsanitized_smiles(stripped_smiles: str) -> Tuple[List, List]:
    """
    Parse unsanitized SMILES into atom and bond lists.
    
    Args:
        stripped_smiles: SMILES string with RDKit atom maps (from Regex stage)
    
    Returns:
        (atoms, bonds)
        - atoms: list of rdkit.Chem.Atom objects
        - bonds: list of (i, j, bond_type) tuples
    
    Raises:
        RDKitException (deferred to Integration phase for re-raise as ValueError)
    """
    
    # Parse with sanitize=False to preserve aromatic flags and skip valence checks
    mol = Chem.MolFromSmiles(stripped_smiles, sanitize=False)
    
    if mol is None:
        raise ValueError(f"Failed to parse unsanitized SMILES: {stripped_smiles}")
    
    # Extract atoms
    atoms = list(mol.GetAtoms())
    
    # Extract bonds
    bonds = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        bonds.append((i, j, bond_type))
    
    # Verify indices are stable (atom indices match GetIdx())
    for atom in atoms:
        assert atom.GetIdx() < len(atoms), f"Atom index {atom.GetIdx()} out of bounds"
    
    # Document implicit hydrogen counts (for verification)
    for atom in atoms:
        implicit_h = atom.GetTotalNumHs(implicitH=True)
        # Log for debugging: f"Atom {atom.GetIdx()} ({atom.GetSymbol()}) has {implicit_h} implicit H"
    
    return atoms, bonds
```

---

**Status:** Ready for Implementation  
**Effort Estimate:** ~30 minutes  
**Blocked By:** MiniPRD_DirectParser_RegexPreprocessor (produces `stripped_smiles` input)
