---
trigger: always_on
glob: "**/*"
description: Standard for Package Management using uv
---

# Package Management Standards (uv)

## 1. The Single Source of Truth
- **Manifest**: `pyproject.toml` is the ONLY source of truth for dependencies.
- **Lockfile**: `uv.lock` must be committed and kept in sync.

## 2. Dependency Management
- **Adding Packages**: ALWAYS use `uv add <package>`.
  - *Bad*: `pip install <package>` (Do not use).
  - *Bad*: Manually editing `pyproject.toml` (unless necessary for config).
- **Dev Dependencies**: Use `uv add --dev <package>` for testing/linting tools.

## 3. Environment
- **Virtual Environment**: `uv` manages the venv specifically for this project.
- **Execution**: Run scripts via `uv run <script>` or ensure the venv is activated.

## 4. Workflows
- **Install**: `uv sync` to install from lockfile.
- **Update**: `uv lock --upgrade` to update dependencies.
