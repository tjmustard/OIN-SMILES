# MiniPRD: Direct Parser — Integration with OIN3DGenerator

**Hypergraph Node ID:** `atom_direct_parser_integration`  
**Parent Node:** `system_oin3d_generator`  
**Blocked By:** All three components: `MiniPRD_DirectParser_RegexPreprocessor`, `MiniPRD_DirectParser_ASTTokenization`, `MiniPRD_DirectParser_MolassemblerInstantiation`  
**Feature:** OIN to Molassembler Direct Parser (v0.2.1)

---

## 1. The Confidence Mandate

**Confidence Score:** 9/10

**Rationale:**
- Three sub-modules (Regex, AST, Molassembler) are independently tested
- Integration is straightforward pipeline: regex → AST → Molassembler → DG
- Timeout handling is TBD (default 10s, hotfix v0.2.2) but doesn't block integration
- Return type (`GeneratedStructure`) already exists and is compatible

**Clarifying Questions:**
1. ✅ **Timeout threshold:** Default 10s hard limit; user feedback informs v0.2.2 adjustment
2. ✅ **DG seeding:** Use default deterministic seed (or 42 if needed); Verification MiniPRD validates
3. ✅ **Error translation:** All exceptions → `ValueError` with context; original preserved in `__cause__`

---

## 2. Atomic User Stories

- **US-001:** As `OIN3DGenerator.generate()`, I want to call the direct parser instead of `_stitch_multi_eta_fragment()`
- **US-002:** As `generate()`, I want to invoke Molassembler's distance geometry engine with 10s timeout
- **US-003:** As `generate()`, I want to return `GeneratedStructure(xyz, mol)` with valid RDKit mol and XYZ block
- **US-004:** As `generate()`, I want to translate low-level Molassembler/RDKit errors into user-friendly `ValueError` messages
- **US-005:** As API, I want to preserve backward compatibility: signature unchanged, return type `GeneratedStructure` per v0.2.0

---

## 3. Implementation Plan (Task List)

- [ ] **Task 3.1:** Replace `_stitch_multi_eta_fragment()` call in `OIN3DGenerator.generate()` with new parser (5 min)
  - Old: `result_xyz, result_mol = _stitch_multi_eta_fragment(oin_string, mol_template)`
  - New: `gen_struct = parse_oin_direct(oin_string)` where `gen_struct.xyz`, `gen_struct.mol`
- [ ] **Task 3.2:** Implement `parse_oin_direct(oin_smiles: str) → GeneratedStructure` wrapper (10 min)
  - Calls: Regex → AST → Molassembler instantiation
  - Error handling: wrap all exceptions as `ValueError`
- [ ] **Task 3.3:** Implement DG conformer generation with timeout: `masm.io.DG(mol, seed=42, timeout_s=10.0)` (10 min)
  - Timeout: use `concurrent.futures.ProcessPoolExecutor` (GIL-safe isolation)
  - Fallback: if timeout, raise `TimeoutError` with message "Conformer generation exceeded 10s; complex geometry may not be supported"
- [ ] **Task 3.4:** Extract XYZ coordinates from Molassembler result; format as XYZ block (5 min)
  - Use `masm.io.write(mol, format='xyz')` or similar
  - Verify XYZ format matches existing output (atomic symbols, coordinates, no metadata)
- [ ] **Task 3.5:** Wrap RDKit mol from AST stage in `GeneratedStructure(xyz, mol)` and return (2 min)
- [ ] **Task 3.6:** Add docstring to `OIN3DGenerator.generate()`: "Internal implementation detail; not part of stable API. Use SMILESToXYZ instead." (5 min)
- [ ] **Task 3.7:** Update error handling: catch all exceptions, re-raise as `ValueError` with clear message (5 min)
- [ ] **Task 3.8:** Integration test: call `generate()` on Cisplatin, TiCat1/3/4, verify output is `GeneratedStructure` (10 min)

---

## 4. The Negative Space (Constraints)

- **DO NOT** change the `OIN3DGenerator.generate()` signature (input/output must remain compatible)
- **DO NOT** log or expose internal Molassembler error messages directly to users (translate to actionable messages)
- **DO NOT** fall back to old `_stitch_multi_eta_fragment()` or `_embed_fragment()` on errors (fail fast)
- **DO NOT** attempt multi-conformer generation or candidate routing (single DG output only)
- **DO NOT** modify the return type from `GeneratedStructure` (v0.2.0 API contract)

---

## 5. Integration Tests & Verification

**Test 1 (Deterministic — Cisplatin through Full Pipeline):**
- **Input:** OIN-SMILES string (v3.6 inline format) `[Pt_SQP].[Cl]{1}.[Cl]{3}`
- **Expected Output:**
  ```python
  gen_struct = OIN3DGenerator.generate(oin_smiles)
  assert isinstance(gen_struct, GeneratedStructure)
  assert gen_struct.xyz is not None  # XYZ block string
  assert gen_struct.mol is not None  # RDKit Mol with conformer
  assert gen_struct.mol.GetNumAtoms() == 5
  ```
- **Verification:**
  - No exception raised
  - XYZ format valid (can be parsed by external tools)
  - Mol has 3D conformer (confirmed by `Chem.Get3DDistanceMatrix()`)

**Test 2 (Deterministic — TiCat1 Ansa-Metallocene):**
- **Input:** Real OIN-SMILES for TiCat1 (ansa-ligand with eta bonds)
- **Expected Output:**
  - `GeneratedStructure(xyz, mol)` where mol has ansa-ring with eta bonds
  - XYZ has Ti-C distances consistent with eta bonding (shorter than vdW sum)
- **Verification:**
  - XYZ generated without error
  - Atom count and bond count match input
  - Round-trip: `xyz2mol(xyz)` → OIN2; compare with original OIN (canonical form)

**Test 3 (Error Case — Timeout on Large Complex):**
- **Input:** Simulated large OIN-SMILES (>100 atoms, complex geometry)
- **Expected Behavior:**
  - DG generation exceeds 10s timeout
  - `TimeoutError("Conformer generation exceeded 10s; complex geometry may not be supported")` raised
  - User can handle with try/except and decide: log warning or fail
- **Verification:**
  ```python
  try:
      gen_struct = generate(large_oin)
  except TimeoutError as e:
      print(f"Warning: {e}")
  ```

**Test 4 (Error Case — Malformed OIN-SMILES):**
- **Input:** `[Pd_UNKNOWN_SHAPE].[Cl]` (invalid shape)
- **Expected Behavior:**
  - Regex extracts `UNKNOWN_SHAPE`
  - Molassembler instantiation fails on shape lookup
  - Exception caught and re-raised as:
    ```python
    ValueError("Failed to construct Molassembler molecule: Unknown shape UNKNOWN_SHAPE; valid: [SQP, OC, TBP, ...]")
    ```
- **Verification:**
  - User gets actionable error message
  - No partial molecule leaked
  - Original exception preserved in `__cause__` for debugging

**Test 5 (Backward Compatibility):**
- **Input:** Same OIN-SMILES as v0.2.0
- **Expected Output:**
  - Return type: `GeneratedStructure` (unchanged)
  - XYZ block: valid and parseable
  - Mol object: consistent with v0.2.0 expectations (may differ in geometry due to DG change, but topology same)
- **Verification:**
  - API contract maintained
  - No breaking changes to callers

---

## Appendix: Code Outline

```python
from concurrent.futures import ProcessPoolExecutor
import scine_molassembler as masm
from typing import Optional

class OIN3DGenerator:
    """
    Generate 3D coordinates for OIN-SMILES strings using Molassembler Distance Geometry.
    
    **Warning:** This class is an internal implementation detail of OIN-SMILES. Its API
    is subject to change without notice. Users should prefer SMILESToXYZ.convert() instead.
    """
    
    @staticmethod
    def generate(oin_smiles: str) -> GeneratedStructure:
        """
        Parse OIN-SMILES and generate 3D coordinates via Molassembler DG.
        
        Args:
            oin_smiles: OIN-SMILES string in v3.6 inline format (e.g., "[Pd_SQP].[Cl]{1}")
        
        Returns:
            GeneratedStructure(xyz, mol) where:
            - xyz: XYZ block string with 3D coordinates
            - mol: RDKit Mol with bond topology and 3D conformer (unsanitized)
        
        Raises:
            ValueError: if OIN parsing fails or Molassembler construction fails
            TimeoutError: if DG conformer generation exceeds 10s timeout
        """
        
        try:
            # Pipeline: Regex → AST → Molassembler → DG
            
            # Step 1: Regex preprocessing
            stripped_smiles, constraints = extract_oin_constraints(oin_smiles)
            
            # Step 2: AST tokenization
            atoms, bonds = tokenize_unsanitized_smiles(stripped_smiles)
            
            # Step 3: Molassembler instantiation
            mol_masm = construct_molassembler_mol(atoms, bonds, constraints, mol_rdkit=None)
            
            # Step 4: DG conformer generation with timeout
            xyz = generate_conformer_with_timeout(mol_masm, timeout_s=10.0)
            
            # Step 5: Package result
            # Note: mol_rdkit from AST stage (unsanitized) is returned alongside XYZ
            mol_rdkit = Chem.MolFromSmiles(stripped_smiles, sanitize=False)  # Rebuild for return
            
            return GeneratedStructure(xyz=xyz, mol=mol_rdkit)
        
        except TimeoutError:
            raise  # Re-raise timeout as-is
        except Exception as e:
            raise ValueError(f"Failed to generate 3D structure from OIN-SMILES: {e}") from e

def generate_conformer_with_timeout(mol_masm: masm.Molecule, timeout_s: float) -> str:
    """
    Generate conformer via Molassembler DG with timeout.
    
    Args:
        mol_masm: Molassembler Molecule
        timeout_s: timeout in seconds (default 10.0)
    
    Returns:
        XYZ block string
    
    Raises:
        TimeoutError: if DG exceeds timeout
    """
    
    def dg_worker(mol):
        # DG conformer generation (runs in separate process)
        try:
            confs = masm.io.generate_conformation(mol, seed=42)
            xyz_str = masm.io.write(confs, format='xyz')
            return xyz_str
        except Exception as e:
            raise RuntimeError(f"DG generation failed: {e}") from e
    
    # Use ProcessPoolExecutor for GIL-safe timeout
    with ProcessPoolExecutor(max_workers=1) as executor:
        future = executor.submit(dg_worker, mol_masm)
        try:
            xyz = future.result(timeout=timeout_s)
            return xyz
        except TimeoutError:
            raise TimeoutError(f"Conformer generation exceeded {timeout_s}s; complex geometry may not be supported")
```

---

**Status:** Ready for Implementation  
**Effort Estimate:** ~55 minutes  
**Blocked By:** MiniPRD_DirectParser_RegexPreprocessor, MiniPRD_DirectParser_ASTTokenization, MiniPRD_DirectParser_MolassemblerInstantiation
