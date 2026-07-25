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

## Handling Hook Rejections
If your `git commit` or `git push` fails:
1. **Commit blocked?** Check if you accidentally staged a file in `spec/process/`. Run `git restore --staged <file>` to fix it.
2. **Push blocked?** Check your branch name. If it starts with `research/` or `swimlane/`, it cannot be pushed. To share the work, you must create a standard feature branch (e.g. `feature/my-change`).
