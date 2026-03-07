# Testing Strategy

## Framework
- **Runner**: `pytest`
- **Manager**: `uv run pytest`

## Test Levels
### 1. Unit Tests (`tests/unit/`)
- **Scope**: Individual functions (`math`, `sorting`, `parsing`).
- **Speed**: Must run in < 100ms.
- **Mocking**: Mock filesystem and heavy RDKit generation if slow.

### 2. Integration Tests (`tests/integration/`)
- **Scope**: Full pipeline (`XYZ -> OIN`, `OIN -> Dict`).
- **Golden Files**: Compare outputs against `tests/data/golden/`.

### 3. Round-Trip Verification (`tests/integration/verify_roundtrip.py`)
- **Mandatory Check**: Must run `bash tests/run_roundtrip.sh` before PR.
- **Criteria**:
  - 100% Success Rate.
  - RMSD < 0.2A (Geometry match).
  - Identity Match (SMILES match).

## Rules
- **No Regression**: Every bug fix requires a test case.
- **Test Data**: Use `.xyz` files in `tests/integration/` (e.g., `Ferrocene.xyz`).
- **CI/CD**: Tests must pass in the `uv` environment.