---
description: Read all memory/rules/docs and entire codebase, then update memory.
---

1.  **Read Baseline**:
    - Read all files in `.agent/memory/`.
    - Read all files in `.agent/rules/`.
    - Read relevant high-level docs in `docs/` (`PRD.md`, `DESIGN.md`).

2.  **Explore Current State**:
    - Use `list_dir` to get a high-level view of the project structure (`src/`, `tests/`, etc.).
    - Use `view_file` to inspect key files (entry points, main logic, config) to understand implementation details.
    - *Goal*: Build a fresh mental model of the *actual* codebase state vs. what is documented.

3.  **Synthesize Findings**:
    - Compare your exploration findings with the current memory state.
    - Identify discrepancies (e.g., deprecated dependencies still mentioned, new patterns not documented).

4.  **Refresh Memory**:
    - **Update `activeContext.md`**: What is the current focus of the codebase? What are the active development areas?
    - **Update `systemPatterns.md`**: Are there new recurring patterns? Have old ones been replaced? Document the *actual* architecture.
    - **Update `productContext.md`**: Are the product goals still aligned with the code? Update if features have evolved.

5.  **Sync Documentation**:
    - If `docs/` are significantly out of date, update `PRD.md` or `DESIGN.md` to match reality.
    - If rules are outdated or new rules are needed, update `.agent/rules/`.

6.  **Verify**:
    - Ensure all memory files are consistent with each other.
    - Ensure memory reflects the codebase accurately.
