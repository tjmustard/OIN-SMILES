# OIN-SMILES - AGENTS.md

## Context
This repository contains the "OIN-SMILES" library.
It facilitates the conversion of 3D transition metal complexes to and from 1D SMILES format.
It follows the MCAF (Managed Code Coding AI Framework).

## Product Requirements
- **Location**: All Product Requirements Documents (PRDs) are located in the `PRDs/` folder.
- **Usage**: Use these documents as the primary source of truth for feature context, technical constraints, and content.
- **Protocol**: Review the relevant PRD in `PRDs/` before drafting implementation plans or modifying code.

## Development Flow
1. **Describe**: Update `docs/Features/` before coding, referencing the specific file in `PRDs/`.
2. **Plan**: Create `implementation_plan.md` and get approval.
3. **Implement**: Write code and tests together.
4. **Verify**: Run tests and static analysis.
5. **Document**: Update docs to reflect reality.

## Maintainer Preferences
- **Language**: Python (managed via `uv`).
- **Style**: Adhere to PEP 8. Use Type Hints for all function signatures.
- **Documentation**: Use clear docstrings for all public modules and functions.
- **Architecture**: Keep parsing logic separated from file I/O.
- **Domain**: Ensure chemical validity in SMILES strings where applicable.

## Testing Discipline
- **Unit Tests**: `uv run python -m unittest discover tests/unit`
- **Integration Tests**: `uv run python -m unittest discover tests/integration`
- **Full Suite**: `uv run python -m unittest discover tests`
- **Linting**: `uv run ruff check .` (or `uv run pylint`, verify preference)

## Commands
- **Install/Sync**: `uv sync`
- **Add Dependency**: `uv add <package>`
- **Test (All)**: `uv run python -m unittest discover tests`
- **Build/Package**: `uv build`

## Self-Learning
- If test discovery fails, verify the `__init__.py` files exist in test directories.
- If import errors occur during `uv run`, ensure `uv sync` has been run to update the virtual environment.
- Document any specific chemical edge cases (e.g., hapticity handling) in `docs/Development/chemistry-edge-cases.md`.