# Design Specification

## Metadata
- **Status**: Computed / Living
- **Scope**: API & CLI Interfaces

## 1. Library Interface (Python API)
### 1.1 Core Canonicalization (`src/oinsmiles/core/translator.py`)
```python
from oinsmiles.core.translator import XYZToSMILES, SMILESToXYZ

# XYZ -> OIN
converter = XYZToSMILES()
oin_string = converter.convert("path/to/file.xyz", charge=0)

# OIN -> Graph (Reconstruction)
reconstructor = SMILESToXYZ()
tmc_graph = reconstructor.convert(oin_string)
```

### 1.2 Generation Engine (`src/oinsmiles/generation/engine.py`)
```python
from oinsmiles.generation.engine import OIN3DGenerator

# OIN -> Architector Dict -> 3D Structure
generator = OIN3DGenerator(scaling_factor=1.0)
structure = generator.generate("[Pt_SPL].[Cl]{0}...", extra_params={"debug": True})
```

## 2. CLI Interface (Proposed)
The system should expose a CLI entry point `oin`.

### 2.1 Commands
- `oin canonicalize <input.xyz> [--charge 0]`
  - Outputs the OIN string to stdout.
- `oin generate <oin_string> [--out output.xyz]`
  - Generates 3D structure using Architector.

## 3. Data Structures
### 3.1 Haptic Vectors
- **Storage**: `TEMPLATES` in `core/geometry_templates.py`.
- **Expansion**: Runtime calculation via `core/haptic_math.expand_slot(n, pos, ref)`.
- **Scaling**: Bond distance heuristic applied in `ArchitectorAdapter`.

## 4. Error Handling
- **Sanitization**: RDKit errors during `MolFromSmiles` must raise `OINSanitizationError`.
- **Geometry**: Failure to match template RMSD raises `GeometryMismatchError`.
