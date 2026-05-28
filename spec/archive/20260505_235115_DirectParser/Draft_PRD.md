# SuperPRD: OIN to Molassembler Direct Parser

## 1. Introduction & Goals

**Problem Statement:**  
The current generation pipeline (`_stitch_multi_eta_fragment`) relies on RDKit's ETKDG for initial 3D generation. This fundamentally fails on ansa-metallocenes because RDKit enforces Hückel's aromaticity rules and cannot kekulize extracted Cp/indenyl rings. Previous workarounds (de-aromatization) corrupted the downstream chemical graph, causing xyz2mol round-trip failures and breaking the OIN stability guarantee.

**Solution Overview:**  
Replace the entire `_stitch_multi_eta_fragment` and `_embed_fragment` logic with a direct parsing architecture that feeds an unsanitized abstract syntax tree (AST) directly into SCINE Molassembler. Molassembler handles haptic bonding natively without enforcing RDKit's aromaticity constraints, preserving chemical graph integrity.

**Target Audience:**  
Internal: The `OIN3DGenerator.generate()` method and downstream 3D conformation pipeline.

---

## 2. Confidence Mandate

**Confidence Score:** 8/10

**Rationale:**  
- Core mutation is well-defined (regex parser → RDKit AST → Molassembler instantiation).
- Molassembler's haptic bond API is proven in production (v0.2.0+).
- Error boundaries are explicit (vertex index validation against polyhedral shapes).
- Verification protocol is deterministic (RMSD < 1.0 Å vs. ground-truth XYZ).

**Open Questions / Clarifications:**  
- None critical. Phase 3 (internal use only) and Phase 4 (deterministic outputs) are resolved.

---

## 3. Scope

**In-Scope:**
- Regex pre-processor to extract polyhedral shape codes, vertex indices, and haptic blocks from OIN-SMILES
- AST tokenization via RDKit (unsanitized SMILES parsing)
- Sequential atom/bond construction in Molassembler
- Stereochemical assignment (SCINE shape injection, haptic bond type assignment)
- Deterministic conformer generation via Molassembler's Distance Geometry engine

**Out-of-Scope:**
- Validation of OIN-SMILES syntax (pre-parser assumes valid format)
- Candidate artifact routing or multi-conformer generation (single deterministic output)
- Support for non-standard polyhedral shapes beyond those in SCINE's shape library
- Modification to OIN format grammar or v3.6 specification

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
|----|-----------|-------------------|----------|
| US-001 | As `OIN3DGenerator`, I want to parse OIN-SMILES strings with eta blocks directly into Molassembler molecules | 1. Regex extracts all polyhedral codes and eta blocks correctly<br>2. Extracted constraints sync to atom index dictionary<br>3. RDKit yields valid AST for unsanitized SMILES | High |
| US-002 | As the parser, I want to instantiate a Molassembler molecule with standard bonds and atoms | 1. Sequential `mol.push_back(atom)` completes without errors<br>2. Standard bonds (single, double, etc.) attach correctly<br>3. Molecule is ready for shape/haptic injection | High |
| US-003 | As the parser, I want to inject polyhedral shapes and haptic bonds into the Molassembler molecule | 1. Metal center receives correct SCINE shape (e.g., SquarePlanar)<br>2. Each eta-listed atom receives eta bond from metal center<br>3. Bond type is BondType.eta as per SCINE spec | High |
| US-004 | As the verifier, I want to validate RMSD of the first coordination sphere against ground-truth XYZ | 1. RMSD ≤ 1.0 Å for all tested complexes (TiCat1/3/4)<br>2. Peripheral/substituent atoms are excluded from RMSD<br>3. Conformer generation completes in reasonable time (< 10s per complex) | High |
| US-005 | As the error handler, I want to reject invalid vertex indices with clear fatal errors | 1. If eta block references index > max_vertices for shape, raise ValueError<br>2. Error message includes shape name, requested index, and max valid index<br>3. No silent fallbacks or heuristic guessing | High |

---

## 5. Technical Specifications

### 5.1 Architecture & Resolved Trade-offs

**Regex Pre-processor:**
- Extract polyhedral shape codes: `_SQP`, `_OC`, `_TBP`, etc. via pattern `_([A-Z0-9]+)`
- Extract vertex indices: `{1}`, `{2,3}` via pattern `\{([0-9,\s]+)\}`
- Extract haptic blocks: `|eta:1,2,3|` via pattern `\|eta:([0-9,\s]+)\|`
- Replace extracted segments with standard RDKit atom maps (`[Pd:1]`) to preserve node tracking
- Store constraints in a Python dictionary keyed by atom index

**Trade-off (Aromaticity Handling):**
- **Rejected:** De-aromatization + RDKit sanitization (corrupts chemical graph)
- **Chosen:** Unsanitized RDKit AST → Molassembler (preserves graph, requires manual shape/bond injection)
- **Rationale:** Molassembler's native haptic bond API is simpler than reverse-engineering RDKit's bond model

**AST Tokenization:**
- Call `Chem.MolFromSmiles(smiles, sanitize=False)` on regex-stripped SMILES
- Extract atoms and bonds directly from the RDKit molecule object
- No aromaticity flags are evaluated or enforced at this stage

**Molassembler Instantiation:**
1. Create empty `scine_molassembler.Molecule()`
2. Iterate atoms, call `mol.push_back(atom)` for each
3. Iterate bonds, call `mol.add_bond(i, j, bond_type)` for standard bonds
4. Query constraint dictionary for metal center's polyhedral shape
5. Assign shape via `mol.set_shape(metal_idx, SCINE_SHAPE)`
6. For each eta-listed atom index, call `mol.add_bond(metal_idx, lig_idx, BondType.eta)`

**Conformer Generation:**
- Call Molassembler's deterministic Distance Geometry engine
- Seeded with fixed seed (or use first result if DG is deterministic by default)
- Return `GeneratedStructure(xyz, mol)` as per v0.2.0+ API contract

### 5.2 System Graph Blast Radius

**Affected Nodes (architecture.yml):**
- `src/oinsmiles/generation/molassembler_adapter.py::_stitch_multi_eta_fragment` — **Fully replaced**
- `src/oinsmiles/generation/molassembler_adapter.py::_embed_fragment` — **Fully replaced**
- `src/oinsmiles/generation/oin_parser.py::OINParser.parse()` — **May integrate regex logic here**
- `src/oinsmiles/generation/engine.py::OIN3DGenerator.generate()` — **Calls new parser (internal interface only)**

**No Changes Required:**
- `src/oinsmiles/utils/xyz2mol.py` (XYZ→OIN pipeline)
- `src/oinsmiles/oin/parser.py` (sidecar-format parser)
- `src/oinsmiles/generation/engine.py` public API (return type already `GeneratedStructure`)

### 5.3 Execution Checklist

- [ ] **MiniPRD_DirectParser_RegexPreprocessor** — Implement regex extractor and constraint dict
- [ ] **MiniPRD_DirectParser_ASTTokenization** — RDKit unsanitized SMILES parsing
- [ ] **MiniPRD_DirectParser_MolassemblerInstantiation** — Atom/bond construction and shape injection
- [ ] **MiniPRD_DirectParser_Integration** — Wire parser into `OIN3DGenerator.generate()`
- [ ] **MiniPRD_DirectParser_Verification** — RMSD validation against TiCat1/3/4 fixtures

### 5.4 API Contracts / Schema

**Input:**
```python
oin_smiles: str  # e.g., "[Pd_SQP].[Cl]{1}.[Cl]{3}.[N]{2}.[N]{4}|eta:1,2,3"
```

**Output:**
```python
GeneratedStructure(
    xyz: str,           # Full XYZ block string
    mol: Chem.Mol       # RDKit molecule with 3D conformer and bond topology
)
```

**Exception:**
```python
ValueError: str  # Raised if eta block references invalid vertex index
# e.g., "Vertex index 5 invalid for SquarePlanar shape (max 3)"
```

### 5.5 Dependencies

- `scine_molassembler >= 2.0.0` (already pinned in production)
- `rdkit >= 2020.09` (already in pyproject.toml)
- `re` (Python standard library)

---

## 6. Negative Constraints

- **DO NOT** silently fall back to de-aromatization or heuristic guessing on invalid vertex indices
- **DO NOT** sanitize the RDKit molecule after AST extraction (preserves haptic topology)
- **DO NOT** modify the OIN format specification or atom map syntax
- **DO NOT** perform candidate conformer routing or multi-conformer generation
- **DO NOT** assume external callers; design for internal `OIN3DGenerator` use only

---

## 7. Risks & Mitigation

| Risk | Mitigation |
|------|-----------|
| **Molassembler API changes** (newer versions incompatible) | Lock `scine_molassembler >= 2.0.0, < 5.0.0` in pyproject.toml; add CI test for version compatibility |
| **Regex extraction fails on malformed OIN-SMILES** | Document that parser assumes pre-validated OIN format; invalid input is a precondition error (not a parser bug) |
| **RMSD exceeds 1.0 Å on novel complexes** | Validation is limited to TiCat1/3/4 fixtures; document scope; add fixture for any new complex requiring tighter tolerance |
| **Conformer generation takes > 10s** | DG is deterministic but may be slow for large ligands; timeout policy is TBD (recommend 30s soft limit, 60s hard limit per complex) |

---

## 8. Success Metrics

1. **Parsing Coverage:** 100% of OIN-SMILES strings with valid polyhedral shapes and vertex indices parse without error
2. **Geometry Fidelity:** RMSD ≤ 1.0 Å for first coordination sphere across all validation fixtures (TiCat1, TiCat3, TiCat4)
3. **Round-Trip Stability:** OIN→XYZ→OIN conversion preserves chemical graph and stereochemistry (zero corruption)
4. **Performance:** Conformer generation completes in < 10 seconds per complex (TBD: measure v0.2.0 baseline)
5. **API Stability:** Internal `OIN3DGenerator` interface remains unchanged; return type is `GeneratedStructure` as per v0.2.0 contract

---

## Appendix: Implementation Notes

**Phase 1 — Regex Pre-processor:**
```python
import re
def extract_oin_constraints(oin_smiles: str) -> tuple[str, dict]:
    # Pattern for polyhedral codes: _SQP, _OC, etc.
    shape_pattern = r'_([A-Z0-9]+)'
    # Pattern for vertex indices: {1}, {2,3}, etc.
    vertex_pattern = r'\{([0-9,\s]+)\}'
    # Pattern for haptic blocks: |eta:1,2,3|
    eta_pattern = r'\|eta:([0-9,\s]+)\|'
    
    # Extract and store; replace with RDKit atom maps
    constraints = {}
    stripped_smiles = oin_smiles
    # ... implementation details in MiniPRD
    return stripped_smiles, constraints
```

**Error Boundary Example:**
```python
def validate_vertex_index(eta_indices: list[int], shape: str):
    max_vertices = SCINE_SHAPE_VERTICES[shape]  # e.g., 4 for SquarePlanar
    for idx in eta_indices:
        if idx >= max_vertices or idx < 0:
            raise ValueError(
                f"Vertex index {idx} invalid for {shape} shape (max {max_vertices - 1})"
            )
```
