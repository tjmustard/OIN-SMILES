---
name: hyper-init
description: Scaffolds a new project by conducting an interview to collect metadata and generating standard repository documentation files from the framework templates. Use when setting up a new repository.
trigger: /hyper-init
---

# Hyper-Init (Repository Scaffolding Agent)

This skill initializes a standard repository structure for a new project within the Hypergraph framework. It interviews the user to collect key variables and then templates out standard repository files (README, CHANGELOG, CONTRIBUTING, etc.).

---

## When to use this skill

- When the user starts a new project and needs to scaffold standard documentation.
- When the user explicitly runs `/hyper-init`.

---

## Step 1 — Interview the User

Use the `AskUserQuestion` tool to collect the following metadata variables. Ask the questions sequentially or in one batch, but do not proceed until you have all of them:

- `{{PROJECT_NAME}}`: The human-readable name of the project (e.g., "Google Docs Annotator").
- `{{REPO_NAME}}`: The machine-readable repository name (e.g., "mcp-gdoc-annotated").
- `{{GITHUB_USERNAME}}`: The GitHub username or organization name.
- `{{AUTHOR_NAME}}`: The author's full name.
- `{{AUTHOR_EMAIL}}`: The author's email address.
- `{{PROJECT_DESCRIPTION}}`: A short, 1-2 sentence summary of what the project does.

---

## Step 2 — Read Templates

Once the variables are collected, read the generic templates from the `.agents/schemas/project-templates/` directory:

- `CHANGELOG.md`
- `CITATION.cff`
- `CODE_OF_CONDUCT.md`
- `CONTRIBUTING.md`
- `DEVELOPMENT.md`
- `README.md`
- `SECURITY.md`
- `pyproject.toml`

---

## Step 2b — Scaffold `pyproject.toml`

Read `.agents/schemas/project-templates/pyproject.toml` and substitute:

- `{{REPO_NAME}}` → the machine-readable repository name from Step 1
- `{{PROJECT_DESCRIPTION}}` → the project description from Step 1
- `{{PACKAGE_NAME}}` → `{{REPO_NAME}}` with hyphens replaced by underscores (e.g., `my-project` → `my_project`)

Write the result to `pyproject.toml` in the project root.

Then run:

```bash
uv add --dev ruff
uv add --dev pytest
```

This installs the required dev tooling and generates `uv.lock`.

---

## Step 3 — Template Replacement and Output

For each of the 7 documentation files (not `pyproject.toml` — already handled above):
1. Replace all placeholder variables (e.g., `{{PROJECT_NAME}}`) with the user's provided answers.
2. Use the Write File tool to save the modified content directly to the project's root directory (e.g., `/README.md`, `/SECURITY.md`).

---

## Step 4 — Verify Output

Run a brief `ls -la` in the root directory to confirm all 8 files have been created successfully (`pyproject.toml` + the 7 documentation files). Inform the user that scaffolding is complete and point them to `/hyper-document` for any future documentation updates.
