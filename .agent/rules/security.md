# Security & Safety

## 1. Input Sanitization (Critical)
- **SMILES Drift Prevention**:
  - ALL incoming SMILES or Mols must pass through `oin/sanitizer.py`.
  - **Zone A Locking**: Metal-binding atoms MUST have `ExplicitHs` set and `NoImplicit` set to True.
- **Untrusted XML/XYZ**: 
  - Validate XYZ format strictly (Atom Count match).
  - Do not execute arbitrary code embedded in input files (pickle, etc.).

## 2. Execution Constraints
- **File System**:
  - Read/Write ONLY within the workspace or explicitly provided output directories.
- **Secrets**:
  - No API keys in source code. Use `.env`.

## 3. Library Safety
- **RDKit**: 
  - Handle `MolFromSmiles` returning `None` gracefully.
  - Wrap complex geometry ops in try/except blocks to prevent crashes on singular matrices.