# SuperPRD: OIN to Molassembler Direct Parser (v0.2.1)

## Metadata
- **Project Name**: OIN-SMILES
- **Feature**: Replace RDKit ETKDG with direct Molassembler Distance Geometry parser
- **Version**: 0.2.1
- **Status**: Resolved — Ready for Implementation
- **Resolution Agent**: hyper-resolve (2026-05-05)

---

## 1. Introduction & Goals

**Problem Statement:**  
The current generation pipeline (`_stitch_multi_eta_fragment`) relies on RDKit's ETKDG for initial 3D generation. This fails on ansa-metallocenes (e.g., TiCat1/3/4) because RDKit enforces Hückel's aromaticity rules and cannot kekulize extracted Cp/indenyl rings. Previous workarounds (de-aromatization) corrupted the downstream chemical graph, causing xyz2mol round-trip failures and breaking the OIN stability guarantee.

**Solution Overview:**  
Replace `_stitch_multi_eta_fragment` and `_embed_fragment` with a direct parsing architecture that feeds an unsanitized RDKit AST directly into SCINE Molassembler. Molassembler handles eta bonding natively without enforcing RDKit's aromaticity constraints, preserving chemical graph integrity.

**Target Audience:**  
**Internal implementation detail.** The `OIN3DGenerator.generate()` method is exported in the public API, but direct use is not recommended. Callers should use `SMILESToXYZ().convert()` instead. **Caveat emptor:** subject to change without notice in future versions.

---

## 2. Confidence Mandate

**Confidence Score:** 8/10

**Rationale:**  
- Core mutation is well-defined (regex parser → RDKit AST → Molassembler instantiation)
- Molassembler's eta bond API is proven in production (v0.2.0+)
- Error boundaries are explicit (vertex index validation against polyhedral shapes)
- Determinism assumed and validated via CI testing
- Trade-offs documented and resolved (timeout, RMSD, validation, versions)

**Red Team Clarifications Resolved:**
1. ✅ Timeout handling → 10s default; hotfix in v0.2.2 post-user-feedback
2. ✅ Determinism guarantee → Assume DG deterministic; validate via CI stability
3. ✅ RMSD baseline → Keep 1.0Å empirical; adjust if test fixtures fail
4. ✅ Input validation → Trap Molassembler errors; re-raise with user-friendly messages
5. ✅ Public API clarity → Documented as internal; kept exported (caveat emptor)
6. ✅ Version compatibility → Pin >=2.0.0, <5.0.0; CI matrix [2.0.0, 3.1.0, 4.8.0]
7. ✅ Unsanitized mol semantics → Documented atomic guarantees; chiral tagging deferred v0.3.0
8. ✅ Error transaction semantics → All-or-nothing mol construction
9. ✅ Candidate artifact routing → Single deterministic DG path; no fallback

---

## 3. Scope

**In-Scope:**
- Regex pre-processor: extract polyhedral shape codes, vertex indices, eta blocks from OIN-SMILES
- AST tokenization via RDKit (unsanitized SMILES parsing)
- Sequential atom/bond construction in Molassembler
- Shape assignment (SCINE shape injection)
- Eta bond type assignment
- Deterministic conformer generation via Molassembler DG engine
- Error handling: trap Molassembler exceptions, re-raise as user-friendly `ValueError`

**Out-of-Scope:**
- Validation of OIN-SMILES syntax (assumes pre-validated; invalid input triggers defensive re-raise)
- Candidate artifact routing or multi-conformer generation (single deterministic output)
- Non-standard polyhedral shapes beyond SCINE library
- Modification to OIN v3.6 format specification
- Chiral tagging of atoms (future work: v0.3.0+)

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
|----|-----------|-------------------|----------|
| US-001 | As `OIN3DGenerator`, I want to parse OIN-SMILES with eta blocks into Molassembler molecules | 1. Regex extracts polyhedral codes and eta blocks correctly<br>2. Constraints sync to atom index dictionary<br>3. RDKit yields valid AST for unsanitized SMILES | High |
| US-002 | As the parser, I want to instantiate Molassembler molecule with atoms and standard bonds | 1. `mol.push_back(atom)` completes without errors<br>2. Standard bonds (single, double, etc.) attach correctly<br>3. Molecule ready for shape/eta injection | High |
| US-003 | As the parser, I want to inject polyhedral shapes and eta bonds | 1. Metal center receives correct SCINE shape (e.g., SquarePlanar)<br>2. Each eta-listed atom receives eta bond from metal<br>3. Bond type is `BondType.eta` per SCINE spec | High |
| US-004 | As verifier, I want to validate RMSD of first coordination sphere | 1. RMSD ≤ 1.0 Å for tested complexes (TiCat1/3/4)<br>2. Peripheral atoms excluded from RMSD<br>3. Conformer generation < 10s per complex (empirical) | High |
| US-005 | As error handler, I want to reject invalid vertex indices with clear errors | 1. Eta block index ≥ max_vertices → raise `ValueError`<br>2. Message includes shape, index, max valid<br>3. No silent fallbacks | High |
| US-006 | As downstream caller, I want robust error on Molassembler failure | 1. Trap low-level exceptions<br>2. Re-raise as `ValueError` with context<br>3. Original exception in `__cause__` | High |

---

## 5. Technical Specifications

### 5.1 Architecture & Resolved Trade-offs

**MVP-First Trade-offs Chosen:**

| Decision | Option | Rationale |
|----------|--------|-----------|
| **Timeout Policy** | Default 10s; hotfix v0.2.2 | Measure baseline post-release; user feedback informs thresholds |
| **Determinism** | Assume; validate in CI | SCINE DG is mature; CI test stability monitors variance |
| **RMSD Baseline** | Keep 1.0Å empirical | Adjust threshold if test fixtures fail; ground-truth measured v0.2.2 |
| **Input Validation** | Trap errors + re-raise | Faster than pre-validation; Molassembler errors guide user |
| **Public API** | Document internal; don't hide | Maintain stability; users opt into risk via docstring warning |
| **Error Atomicity** | All-or-nothing mol | No partial results; prevents silent corruption |
| **Fallback Strategy** | Single DG path; no template fallback | Removes complexity; determinism trumps coverage |

**Regex Pre-processor:**
```
Shape codes:   _([A-Z0-9]+)           e.g., _SQP, _OC
Vertex indices: \{([0-9,\s]+)\}       e.g., {1}, {2,3}
Eta blocks:    \|eta:([0-9,\s]+)\|    e.g., |eta:1,2,3|
```
- Extract and store in constraint dict (keyed by atom index)
- Replace extracted segments with RDKit atom maps ([Pd:1]) for tracking
- No validation of format (assumes valid OIN input)

**AST Tokenization:**
- Call `Chem.MolFromSmiles(smiles, sanitize=False)` on regex-stripped SMILES
- Extract atoms and bonds directly from RDKit mol
- **Guarantee:** Implicit hydrogens preserved; valence not checked; aromaticity not enforced

**Molassembler Instantiation (All-or-Nothing):**
1. Create empty `scine_molassembler.Molecule()`
2. Iterate atoms → `mol.push_back(atom)` for each
3. Iterate bonds → `mol.add_bond(i, j, bond_type)` for standard bonds
4. Query constraints for metal shape
5. `mol.set_shape(metal_idx, SCINE_SHAPE)`
6. For each eta index → `mol.add_bond(metal_idx, lig_idx, BondType.eta)`
7. **Error handling:** Any failure → catch, re-raise as `ValueError` with context; don't return partial mol

**Conformer Generation:**
- Call Molassembler DG engine
- Seed: 42 (or let DG use internal seed if deterministic by default)
- Timeout: **10 seconds hard limit** (may increase v0.2.2)
- Return `GeneratedStructure(xyz, mol)` per v0.2.0 API

### 5.2 System Graph Blast Radius

**Affected Nodes:**
- `molassembler_adapter.py::_stitch_multi_eta_fragment` — **Replaced**
- `molassembler_adapter.py::_embed_fragment` — **Replaced**
- `oin_parser.py::OINParser.parse()` — **May integrate regex logic**
- `engine.py::OIN3DGenerator.generate()` — **Calls new parser**

**No Changes:**
- `perception_tmc.py` (XYZ→OIN pipeline)
- `oin/parser.py` (sidecar format)
- `engine.py` public API (return type already `GeneratedStructure`)

### 5.3 Execution Checklist (MiniPRDs)

- [ ] **MiniPRD_DirectParser_RegexPreprocessor** — Regex extractor + constraint dict
- [ ] **MiniPRD_DirectParser_ASTTokenization** — RDKit unsanitized parsing
- [ ] **MiniPRD_DirectParser_MolassemblerInstantiation** — Atom/bond/shape construction
- [ ] **MiniPRD_DirectParser_Integration** — Wire into `OIN3DGenerator.generate()`
- [ ] **MiniPRD_DirectParser_Verification** — RMSD validation (TiCat1/3/4)

### 5.4 API Contracts

**Input:**
```python
oin_smiles: str
# e.g., "[Pd_SQP].[Cl]{1}.[Cl]{3}.[N]{2}.[N]{4}|eta:1,2,3"
```

**Output (Success):**
```python
GeneratedStructure(
    xyz: str,           # Full XYZ block with 3D coordinates
    mol: Chem.Mol       # RDKit mol with bonds and conformer
)
```

**Unsanitized Mol Guarantees:**
- ✅ Atom count, type preserved from OIN-SMILES
- ✅ Bonds (single, double, aromatic, eta)
- ✅ 3D conformer from Molassembler DG
- ❌ **NOT** kekulized (aromatic flags may differ)
- ❌ **NOT** valence-checked (implicit H not exact)
- ❌ **NOT** chiral-tagged (no @/@@ in `.GetSmiles()`) — *Future: v0.3.0+*

**Exception:**
```python
ValueError(str)
# "Failed to construct Molassembler molecule: Vertex index 5 invalid for SquarePlanar (max 3)"
# Original exception in __cause__
```

### 5.5 Dependencies & Version Compatibility

**Pinned:**
- `scine_molassembler >= 2.0.0, < 5.0.0`
- `rdkit >= 2020.09`
- `re` (stdlib)

**CI Test Matrix:** [2.0.0, 3.1.0, 4.8.0]  
**Known Limitation:** Versions >4.8.0 untested; API drift possible

---

## 6. Negative Constraints

- **DO NOT** silently fall back on errors
- **DO NOT** sanitize RDKit mol after AST extraction
- **DO NOT** modify OIN v3.6 specification
- **DO NOT** perform multi-conformer routing
- **DO NOT** return half-constructed molecules

---

## 7. Risks & Mitigation

| Risk | Mitigation |
|------|-----------|
| **Molassembler API drift** | CI matrix test on [2.0.0, 3.1.0, 4.8.0]; document tested versions |
| **DG timeout > 10s** | Default 10s; measure baseline v0.2.2; increase if needed |
| **DG variance across hardware** | Assume deterministic; monitor CI for flakes |
| **RMSD > 1.0Å on novel complexes** | Validation limited to TiCat1/3/4; add fixture if needed |
| **Unsanitized mol breaks downstream** | Document guarantees; callers must wrap in `Chem.SanitizeMol()` |
| **Malformed OIN crashes** | Trap exceptions; re-raise with user-friendly message |

---

## 8. Success Metrics

1. **Parsing Coverage:** 100% of OIN-SMILES from `XYZToSMILES.convert()` v0.2.0 fixtures parse without error
2. **Geometry Fidelity:** RMSD ≤ 1.0 Å for TiCat1/3/4 first coordination sphere
3. **Round-Trip:** OIN→XYZ→OIN preserves chemical graph (canonical OIN match)
4. **Performance:** Conformer generation < 10s per complex
5. **Error Clarity:** Malformed OIN raises `ValueError` with actionable message
6. **API Stability:** `OIN3DGenerator.generate()` signature unchanged; return `GeneratedStructure`

---

## 9. Future Work (v0.3.0+)

- Chiral tagging: Add @/@@ to atoms via CIP
- Pre-validation: `validate_oin_syntax()` for fail-fast
- Timeout tuning: Measure DG baseline; set conservative thresholds
- Multi-conformer: Candidate artifact routing if needed
- Unsanitized wrapping: Auto-sanitization on return

---

**Resolution Date:** 2026-05-05  
**Status:** Ready for MiniPRD implementation
