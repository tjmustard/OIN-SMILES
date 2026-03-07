---
description: Update agentic memory (docs, rules, memory) with work done in the current session.
---

# /session-update Workflow

1.  **Review Session Work**:
    - Check `task.md` to see what was accomplished.
    - Check git status or recall the file changes made during this session.
    - Identify any key decisions, architectural changes, or new patterns that emerged.

2.  **Update Product & Design Specs**:
    - If Requirements changed: Update `spec/compiled/SuperPRD.md`.
    - If Design changed: Update `docs/DESIGN.md`.
    - *Goal*: Keep high-level docs in sync with the code reality. **CRITICAL:** Do NOT manually modify `spec/compiled/architecture.yml`.

3.  **Update Agentic Memory (Heuristics)**:
    - **Active Context**: Update `.agents/memory/activeContext.md` with the current focus, recent changes, and next steps.
    - **System Patterns**: Update `.agents/memory/systemPatterns.md` if you discovered or created new code patterns, idioms, or architectural structures.
    - **Product Context**: Update `.agents/memory/productContext.md` if there are new product-level decisions.
    - *Note:* These files serve as heuristic context to aid the Builder Agent, acting alongside the deterministic constraints of the Hypergraph.

4.  **Update Rules (Optional)**:
    - If you found yourself repeating a correction or following a new rule, add it to `.agents/rules/active.md` (or creating a new rule file in `.agents/rules/` if it's a major topic).

5.  **Commitment**:
    - Commit the changes to memory and docs alongside the code changes, or as a separate "chore: update memory" commit.
