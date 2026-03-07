# System Patterns
## Purpose
To document the architectural decisions and standards of OIN-SMILES.

## Architecture
**Transformation Bridge Pattern**
- **Ingestion**: `XYZ` -> `RDKit Mol` (via `xyz2mol`).
- **Canonicalization**: `Mol` -> `OIN String` (via `OIN_LIB`).
- **Reconstruction**: `OIN String` -> `Architector inputDict` (via `OIN_LIB`).
- **Generation**: `inputDict` -> `Architector` -> `XYZ` (via External Agent).

## Tech Stack
- **Language**: Python 3.10+
- **Dependency Manager**: `uv`
- **Core Libraries**: `rdkit`, `numpy`, `scipy`
- **External**: `xyz2mol` (vendored), `Architector` (runtime dependency).

## Design Patterns
- **Adapter**: `ArchitectorAdapter` converts internal `ParsedOIN` to external `inputDict` schema.
- **Facade**: `OIN3DGenerator` provides a simple interface for the generation subsystem.
- **Template Method**: `TEMPLATES` logic for geometry matching.
- **Strategy**: Sort strategies (Mass-First Waterfall) for canonicalization.

## Data Contracts
- **OIN String**: `[Metal_GEO].L1{Tags}.L2{Tags}`
- **Architector Dict**: Strict schema with "user_core" and expanded haptic `coordList`.
