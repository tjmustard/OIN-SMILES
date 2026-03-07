---
description: Update agentic memory (docs, rules, memory) with work done in the current session.
---

1.  **Review Session Work**:
    - Check `task.md` to see what was accomplished.
    - Check git status or recall the file changes made during this session.
    - Identify any key decisions, architectural changes, or new patterns that emerged.

2.  **Update Product & Design Docs**:
    - If Requirements changed: Update `docs/PRD.md`.
    - If Design changed: Update `docs/DESIGN.md`.
    - If Plan changed: Update `docs/plan.md`.
    - *Goal*: Keep high-level docs in sync with the code reality.

3.  **Update Agentic Memory**:
    - **Active Context**: Update `.agent/memory/activeContext.md` with the current focus, recent changes, and next steps.
    - **System Patterns**: Update `.agent/memory/systemPatterns.md` if you discovered or created new code patterns, idioms, or architectural structures.
    - **Product Context**: Update `.agent/memory/productContext.md` if there are new product-level decisions or constraints.

4.  **Update Rules (Optional)**:
    - If you found yourself repeating a correction or following a new rule, add it to `.agent/rules/active.md` (or creating a new rule file in `.agent/rules/` if it's a major topic).

5.  **Commitment**:
    - Commit the changes to memory and docs alongside the code changes, or as a separate "chore: update memory" commit.
