
# Product Requirements Document (PRD)

## Metadata
- **Project Name**: OIN-SMILES (Open Isomer Notation & Engine)
- **Version**: 4.0.0 (Master Build)
- **Status**: Approved / Implementation Ready
- **Owner**: Antigravity (derived from User Master Specification)

## 1. Introduction & Goals
### 1.1 Problem Statement
Computational organometallic chemistry lacks a standard, reversible 1D representation for complex 3D stereochemistry. Standard SMILES strings fail to distinguish between stereoisomers (cis/trans/fac/mer) or accurately define haptic bonding orientations (e.g., Ferrocene).

### 1.2 Solution Overview
A **Hashing & Serialization Engine** that creates a **1:1 map between 3D structures and human-readable OIN strings**.
- **Input**: Raw 3D XYZ coordinates.
- **Process**: Canonicalization (sorting, geometry fitting, haptic/chelate tagging).
- **Output**: Unique OIN String (e.g., `[Pt_SPL].[Cl]{0}.[Cl]{1}`).
- **Reverse**: Reconstruction of 3D configuration dictionaries for **Architector**.

### 1.3 Target Audience
- **The Pipeline Engineer**: Needs to canonicalize millions of XYZ structures for deduplication.
- **The Generative Chemist**: Needs to write simple strings to auto-generate complex 3D inputs.

## 2. Confidence Mandate
**Confidence Score**: 10/10 (Master Specification provided)
**Clarifying Questions**:
- [x] Are we architecting TMC-Builder or OIN? (User confirmed OIN-SMILES v4.0).

## 3. Scope
### 3.1 In-Scope
- **Canonicalization Engine**: `xyz2mol` -> Geometry Detection -> OIN String.
- **Reconstruction Engine**: OIN String -> `inputDict` -> Architector.
- **Haptic Logic**: "North Star" vector logic for contiguous rings (Cp, Allyl).
- **Chelate Logic**: Explicit integer mapping `{n}..{m}` for bidentate/tridentate ligands.
- **Sanitization**: RDKit Zone A locking to prevent "SMILES Drift".

### 3.2 Out-of-Scope
- **New Quantum Mechanics Methods**: We use existing methods (GFN2-xTB) via Architector; we do not invent new ones.
- **GUI Development**: CLI and Python API only.

## 4. User Stories (Atomic)
| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As a Pipeline Engineer, I want to convert an XYZ file to a unique OIN string, so that I can deduplicate my database. | 1. Same 3D structure (rotated) -> Same String.<br>2. Different stereoisomer -> Different String. | High |
| US-002 | As a Chemist, I want to write a string like `[Pt_SPL].[Cl]{0}..`, so that I can generate the 3D structure without manual coordinate entry. | 1. String parses correctly.<br>2. Generates valid Architector `inputDict`. | High |
| US-003 | As a Developer, I need the system to handle Ferrocene (Cp2Fe), so that haptic ligands are accurately represented. | 1. Detects `[cH]` rings.<br>2. Assigns `{n>}` or `{n<}` tags based on winding. | High |
| US-004 | As a user, I want mass-based sorting of ligands, so that the string is deterministic. | 1. Heavier fragments listed first.<br>2. Heavier binding atoms listed second. | Medium |

## 5. Technical Specifications (The Blueprint)
### 5.1 Architecture
**Transformation Bridge**:
1.  **Ingestion**: `XYZ` -> `xyz2mol` -> `RDKit Mol`.
2.  **Canonicalization**: `OIN_LIB` (Sanitizer -> Sort -> Geometry Fit -> North Star) -> `OIN String`.
3.  **Reconstruction**: `OIN String` -> `OIN_LIB` (Parse -> Haptic Expand) -> `Architector inputDict`.
4.  **Generation**: `inputDict` -> `Architector` -> `New XYZ`.

### 5.2 API Contracts / Schema
#### The OIN String Syntax
`[Metal_GEO].L1{Tags}.L2{Tags}`
- **GEO**: LIN, SPL, TET, OCT, TBP, TPY.
- **Haptic Tags**: `{n>}` (Forward), `{n<}` (Reverse).
- **Chelate Tags**: `{n}..{m}`.

#### Architector `inputDict` Schema
```python
{
    "core": {
        "metal": "Pt",
        "coreType": "user_core",
        "coordList": [[x,y,z], ...] # Expanded haptic vectors
    },
    "ligands": [
        {
            "smiles": "[Cl]",
            "coordinating_atoms": [0],
            "coordList": [[0, 0]] # Map atom 0 to core vector 0
        }
    ],
    "parameters": { ... }
}
```

### 5.3 Dependencies
- **Core**: `numpy`, `scipy`, `rdkit`.
- **Management**: `uv`.
- **External**: `Architector` (assumed downstream), `xyz2mol` (local/vendored).

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** allow "SMILES Drift" (e.g., `[cH]` optimizing to `c`). **MUST** use Zone A locking.
- **DO NOT** use implicit hydrogens for Zone A atoms.
- **DO NOT** rely on input atom order; **MUST** use Mass-First Waterfall sort.
- **DO NOT** map to the *nearest* slot without first locking the global frame via RMSD.

## 7. Risks & Mitigation
- **Risk 1**: `xyz2mol` graph interpretation failure. -> **Mitigation**: Users must provide valid, non-broken XYZ geometries.
- **Risk 2**: Architector generation failure. -> **Mitigation**: We only guarantee valid `inputDict` generation, not successful QM convergence.

## 8. Success Metrics
- **Round-Trip Fidelity**: 100% of valid inputs recover their stereochemistry after XYZ->OIN->XYZ cycle.
- **Uniqueness**: 0% collision rate for distinct isomers.
