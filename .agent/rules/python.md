---
trigger: always_on
glob: "**/*.py"
description: Python standards for OIN-SMILES project
---

# Python Scientific Standards (OIN-SMILES)

## Tech Stack
- **Language**: Python 3.10+
- **Manager**: `uv` (Strict).
- **Linter**: `ruff` (Strict).
- **Core Libs**: `numpy`, `scipy`, `rdkit`.

## Code Style
- **Type Hints**: REQUIRED for ALL function arguments and return values.
  - `def func(x: List[float]) -> np.ndarray:`
- **Docstrings**: Google Style required.
- **Paths**: `pathlib.Path` ONLY. No string concatenation for paths.

## Architectural Constraints
- **Transformation Bridge**:
  - **Canonicalization**: Pure functions. No side effects. Deterministic.
  - **Generation**: Strict schema compliance for `inputDict`.
- **Dependency Isolation**:
  - `core/` should NOT import from `generation/`.
  - `utils/` should remain generic.

## Anti-Patterns
- **Global State**: Mutable global variables are FORBIDDEN.
- **Implicit Conversions**: Do not rely on implicit RDKit modifications. Use explicit Sanitizers.
- **Print Debugging**: Use `logging` or return structured errors.
- **Bare Exceptions**: Never use `except:`. Catch specific errors (`ValueError`, `ImportError`).