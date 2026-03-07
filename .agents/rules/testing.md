# Testing Strategy

## Framework
- **Runner**: `pytest`
- **Manager**: `uv run pytest`

## Test Levels
### 1. Unit Tests (`tests/unit/`)
- **Scope**: Individual functions and classes.
- **Speed**: Must run fast.
- **Mocking**: Mock external dependencies, filesystem, and heavy computations.

### 2. Integration Tests (`tests/integration/`)
- **Scope**: Full pipeline or multi-component interactions.
- **Golden Files**: Compare outputs against expected data in `tests/data/`.

### 3. Verification
- **Mandatory Check**: Must run verification suite before PR.
- **Criteria**:
  - 100% Success Rate.

## Rules
- **No Regression**: Every bug fix requires a test case.
- **CI/CD**: Tests must pass in the `uv` environment.
