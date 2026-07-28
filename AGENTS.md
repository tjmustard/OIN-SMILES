> **HACF as a Toolchain:** This project uses the Hypergraph Coding Agent Framework
> (HACF) as its development toolchain. The skills in `.agents/skills/`, the scripts
> in `.agents/scripts/`, and the schemas in `.agents/schemas/` are development tools —
> they are **not** subjects of this project's plans, PRDs, or architecture docs.
> When you create SuperPRDs, MiniPRDs, or architecture nodes, you are documenting
> **this project**, not the HACF framework itself.

# OIN-SMILES - AGENTS.md

## Context
This repository contains the "OIN-SMILES" library.
It facilitates the conversion of 3D transition metal complexes to and from 1D SMILES format.
It uses the Hypergraph Coding Agent Framework (HACF) as its development toolchain.

## Product Requirements
- **Location**: All Product Requirements Documents (PRDs) are located in the `PRDs/` folder.
- **Usage**: Use these documents as the primary source of truth for feature context, technical constraints, and content.
- **Protocol**: Review the relevant PRD in `PRDs/` before drafting implementation plans or modifying code.

## Development Flow
1. **Describe**: Write the feature description before coding, referencing the specific file in `PRDs/`.
2. **Plan**: Create `implementation_plan.md` and get approval.
3. **Implement**: Write code and tests together.
4. **Verify**: Run tests and static analysis.
5. **Document**: Update docs to reflect reality — see **Documentation Layout** below for *where*.

## Documentation Layout
`docs/` is split by audience, and **the root is closed to new files**:
- **Product documentation** — `docs/` root. Four files only: `README.md`, `OPTIMIZERS.md`,
  `GENERATION_PIPELINE.md`, `KNOWN_LIMITATIONS.md`. Adding a fifth needs maintainer sign-off.
- **Agentic coding notes** — `docs/agentic-notes/<release>/`. Everything a session produces:
  measurement reports, lane write-ups, A/B results, refuted hypotheses, status snapshots.

**The test:** if it records what a session *measured, tried, or refuted*, it is a note. If it
tells someone *how the shipped software behaves*, it is a product doc. `<release>` is the release
the work is **for** — targeting v0.4.7 means creating `docs/agentic-notes/v0.4.7/`, not reusing
the newest folder that already exists.

- **Frozen measurement data** — `measurements/<release>/`. The numbers a *later release diffs*:
  frozen baselines, gate/audit tallies, transition matrices. Written by
  `tools/harvest_measurements.py`, never by hand. **The test:** *"will a future release diff this
  to know whether it regressed?"* — if yes it is a measurement; if it is the evidence for one
  paragraph of one write-up, it is a note and belongs beside that note.

A `pre-commit` guard (`tools/check_docs_layout.sh`) rejects new files at the `docs/` root. Full
rule: `.agents/rules/docs-layout.md`.

## Maintainer Preferences
- **Language**: Python (managed via `uv`).
- **Style**: Adhere to PEP 8. Use Type Hints for all function signatures.
- **Documentation**: Use clear docstrings for all public modules and functions.
- **Architecture**: Keep parsing logic separated from file I/O.
- **Domain**: Ensure chemical validity in SMILES strings where applicable.

## Version Control
- **Landing**: feature work lands on `main` by **squash-merge** (PR #2/#3 precedent). Never hand-merge a branch's granular history.
- **Tag before you delete.** A squash-merge creates a *new* commit, so the branch's own commits are never ancestors of `main` — `git branch -d` will refuse and `git branch -D` throws the history away. Before deleting, tag the tip:
  ```bash
  git tag -a archive/<slug> -m "Granular history of <session>; content landed as <squash-sha>" <branch>
  git branch -D <branch>
  ```
  Without the tag the commits are left dangling and `git gc` prunes them. Both S3 and S4 lost their history this way; S4's was recovered only because the dangling commit had not yet been collected (see `archive/s3-aromatic-perception`, `archive/s4-eta-winding`).
- **"N commits not in main" after a squash is expected, not a lost merge.** Verify content, not commit ancestry, before deleting: `git diff --quiet main..<branch> -- <path>` per delivered file. Files identical on both sides do **not** appear in `git diff --stat`, so a large-looking stat almost always means `main` is *newer*, not that content is missing.
- **`main` is intentionally ahead of `origin/main` and unpushed.** Do not push or open a PR without asking the maintainer.
- **`git stash push <path>` on an already-committed (clean) path saves nothing and exits 0.** The following `git stash pop` then applies whatever stash *is* on top. To A/B a file against the previous commit use `git show HEAD~1:<path> > <path>` … `git checkout -- <path>`.

## Branching Safeguards & Internal Files
Review `.agents/rules/git-workflow.md` for strict rules regarding protected paths (e.g. `spec/process/`) and blocked remote branches (`research/*`, `swimlane/*`). Agents must respect these constraints.

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
- Document any specific chemical edge cases (e.g., hapticity handling) in `docs/KNOWN_LIMITATIONS.md` — that is the product doc for reproducible failure modes. The *investigation* behind one goes in `docs/agentic-notes/<release>/`.