---
name: session-update
description: Updates agentic memory files (activeContext, systemPatterns, productContext) and project docs with work done in the current session. Use at the end of a development session to keep memory in sync with the codebase.
---

# Session Update

This skill closes out a development session by syncing memory files, documentation, and heuristic context with the work completed. It ensures the next agent session starts with accurate context.

## When to use this skill

- At the end of every development session, before closing the chat.
- When the user explicitly runs `/hyper-session-update`.
- After completing a significant feature, fix, or architectural change.

## How to use it

1. **Review Session Work**
   - Check task tracking docs (e.g., `task.md`) to identify what was accomplished.
   - Review git status or recall the file changes made during this session.
   - Identify key decisions, architectural changes, or new patterns that emerged.

2. **Update Product & Design Specs**
   - If requirements changed: Update `spec/compiled/SuperPRD.md`.
   - If design changed: Update `docs/DESIGN.md`.
   - *Goal:* Keep high-level docs in sync with code reality.
   - **CRITICAL:** Do NOT manually modify `spec/compiled/architecture.yml`. Use the deterministic scripts or `/hyper-discover`/`/hyper-audit` for that.

3. **Update Agentic Memory**
   - **`activeContext.md`**: Update with current focus, recent changes, and immediate next steps.
   - **`systemPatterns.md`**: Add new code patterns, idioms, or architectural structures discovered. Remove replaced ones.
   - **`productContext.md`**: Update if there are new product-level decisions or evolved features.
   - *Note:* These files act as heuristic context alongside the deterministic Hypergraph.

4. **Update Rules (if needed)**
   - If a correction was repeated or a new rule consistently applied, add it to `.agents/rules/` (or create a new rule file for major topics).

5. **Commit**
   - Commit memory and doc updates alongside code changes, or as a separate `chore: update memory` commit.
