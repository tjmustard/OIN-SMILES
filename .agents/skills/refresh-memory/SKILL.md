---
name: refresh-memory
description: Reads all memory, rules, docs, and the codebase to rebuild a fresh mental model of the project, then updates the memory files to reflect actual current state. Use when starting a new session or when memory feels stale.
---

# Refresh Memory

This skill rebuilds the agent's mental model of the project by reading deterministic state, exploring the actual codebase, and updating the heuristic memory files to reflect current reality.

## When to use this skill

- When starting a new session after a significant gap in development.
- When the agent's responses seem to reference outdated architecture or patterns.
- When the user explicitly runs `/refresh-memory`.
- After major refactors that may have drifted from documented patterns.

## How to use it

1. **Read Baseline State**
   - Read the deterministic state: `spec/compiled/architecture.yml` (the Hypergraph).
   - Read the compiled specifications: `spec/compiled/SuperPRD.md` and all `MiniPRD_*.md` files.
   - Read all files in `.agents/memory/` and `.agents/rules/`.

2. **Explore the Current Codebase**
   - Use directory listing to get a high-level view of the project structure (`src/`, `tests/`, etc.).
   - Use search and file reading to inspect key files and understand implementation details.
   - *Goal:* Build a fresh mental model of the *actual* codebase state vs. what the Hypergraph currently documents.

3. **Synthesize Findings**
   - Compare your exploration findings with the current memory state.
   - Identify discrepancies (e.g., deprecated dependencies still mentioned, new patterns not documented, renamed files).

4. **Refresh the Memory Files**
   - **`activeContext.md`**: Update with the current development focus, recent changes, and immediate next steps.
   - **`systemPatterns.md`**: Update with new recurring patterns. Remove replaced ones. Document the *actual* architecture.
   - **`productContext.md`**: Update if product goals have evolved or features have changed.

5. **Sync Documentation**
   - If abstract patterns or rules are outdated, update `.agents/memory/systemPatterns.md` or the relevant `.agents/rules/` file.
   - **CRITICAL:** Do NOT manually edit `architecture.yml`. If the Hypergraph is out of sync with reality, instruct the user to run `/discover` or `/audit` to deterministically reconcile it.

6. **Verify Consistency**
   - Ensure all memory files are consistent with each other and accurately reflect the codebase.
   - Report any discrepancies found that could not be resolved without running a deterministic script.
