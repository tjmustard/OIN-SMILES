# Git & Workflow Standards (OIN-SMILES)

## Protected Paths (DO NOT COMMIT)
The following directories are for internal AI framework use only. They are ignored by git and guarded by a `pre-commit` hook. **Do not attempt to commit files in these directories:**
- `spec/process/`
- `spec/handoffs/`
- `scratchpad/`

## Branching Strategy & Push Protections
- **`main`**: The stable production branch. Feature work lands here via squash-merge.
- **Experimental Branches (`research/*`, `swimlane/*`)**: These prefixes are heavily guarded. A `pre-push` hook physically prevents these branches from ever being pushed to the remote GitHub repository. 
  - **Usage:** Agents MUST use these prefixes for isolated work, refactoring, or running scripts. 
  - **Why:** This ensures experimental or broken agent code never pollutes the remote repository.

## Commit Trailers & `--no-verify`
- A `commit-msg` hook standardizes the agentic-tool session trailer: it rewrites any
  `Claude-Session: <url>` (and the Cursor/Antigravity equivalents) line into a single
  `Agentic Coding Assisted: Claude Code` trailer. **Do not hand-write that trailer** — write the
  session line as normal and let the hook normalize it.
- **Do NOT use `git commit --no-verify` for routine commits.** `--no-verify` bypasses *every*
  client-side hook, including `commit-msg`, so the raw session URL lands in history un-rewritten
  (and it also skips the lint check). The `pre-commit` lint runs `uvx ruff@0.15.20` (isolated,
  ephemeral) — it does **not** re-sync the project env or float the pinned `rdkit==2025.9.3` — so
  there is no longer any reason to bypass hooks for a normal commit.
- If you ever must `--no-verify` deliberately (e.g. a pathspec/collision-safe land in a shared
  worktree), amend afterward (or hand-add the trailer) so history stays standardized.

## Handling Hook Rejections
If your `git commit` or `git push` fails:
1. **Commit blocked?** Check if you accidentally staged a file in `spec/process/`. Run `git restore --staged <file>` to fix it.
2. **Push blocked?** Check your branch name. If it starts with `research/` or `swimlane/`, it cannot be pushed. To share the work, you must create a standard feature branch (e.g. `feature/my-change`).
