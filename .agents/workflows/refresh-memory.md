---
description: Read all memory/rules/docs and entire codebase, then update memory.
---
# /refresh-memory Workflow

1.  **Read Baseline**:
    - Read the deterministic state: `spec/compiled/architecture.yml` (The Hypergraph).
    - Read the active specifications: `spec/compiled/SuperPRD.md` and `MiniPRD`s.
    - Read heuristic memory: all files in `.agents/memory/` and `.agents/rules/`.

2.  **Explore Current State**:
    - Use `list_dir` to get a high-level view of the project structure (`src/`, `tests/`, etc.).
    - Use `grep_search` or `view_file` to inspect key files to understand implementation details.
    - *Goal*: Build a fresh mental model of the *actual* codebase state vs. what the Hypergraph currently documents.

3.  **Synthesize Findings**:
    - Compare your exploration findings with the current memory state.
    - Identify discrepancies (e.g., deprecated dependencies still mentioned, new patterns not documented).

4.  **Refresh Memory**:
    - **Update `activeContext.md`**: What is the current focus of the codebase? What are the active development areas?
    - **Update `systemPatterns.md`**: Are there new recurring patterns? Have old ones been replaced? Document the *actual* architecture.
    - **Update `productContext.md`**: Are the product goals still aligned with the code? Update if features have evolved.

5.  **Sync Documentation**:
    - If abstract patterns or rules are outdated, update `.agents/memory/systemPatterns.md` or `.agents/rules/`.
    - **CRITICAL:** Do NOT manually edit `architecture.yml`. If the Hypergraph is out of sync with reality, instruct the user to run `/discover` or `/audit` to deterministically reconcile it.

6.  **Verify**:
    - Ensure all memory files are consistent with each other.
    - Ensure memory reflects the codebase accurately.
