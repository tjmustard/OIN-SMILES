---
name: document
description: Updates project documentation (README, CHANGELOG, design docs) after code changes. Reads actual code before writing — never trusts existing docs. Use after completing a feature, fix, or refactor.
---

# Document

This skill updates project documentation after code changes. It reads the actual implementation before writing anything — existing documentation is treated as potentially stale.

## When to use this skill

- After completing a feature, bug fix, or refactor that changed user-facing behavior.
- When the user explicitly runs `/document`.
- When README, CHANGELOG, or design docs are out of sync with the codebase.

## How to use it

1. **Identify Changes**
   - Check `git diff` or recent commits for modified files.
   - Identify which features/modules were changed.
   - Note any new files, deleted files, or renamed files.

2. **Verify Current Implementation**
   **CRITICAL:** Do NOT trust existing documentation. Read the actual code.
   - Read the current implementation of each changed file.
   - Understand actual behavior, not documented behavior.
   - Note any discrepancies between existing docs and actual behavior.

3. **Update Relevant Documentation**

   - **`README.md`**: Reflect architectural, structural, or high-level project workflow changes.
   - **`CHANGELOG.md`**:
     - **CRITICAL:** Before modifying the CHANGELOG, ask the user: *"Are we bumping the version? If so, what is the new semantic version and today's date?"*
     - Wait for their response.
     - **If a version/date is provided:** Rename `## [Unreleased]` to `## [Version] - Date` and create a fresh `## [Unreleased]` block above it.
     - **If no version bump:** Add entries under the existing `## [Unreleased]` section.
     - Use categories: `Added`, `Changed`, `Fixed`, `Security`, `Removed`.
     - Write in concise, user-facing language.

4. **Documentation Style Rules**

   ✅ **Concise** — Sacrifice grammar for brevity
   ✅ **Practical** — Examples over theory
   ✅ **Accurate** — Code-verified, not assumed
   ✅ **Current** — Matches actual implementation

   ❌ No enterprise fluff
   ❌ No outdated information
   ❌ No assumptions without verification

5. **Ask if Uncertain**
   If the intent behind a change or its user-facing impact is unclear, ask the user — don't guess.
