---
description: Update documentation after code changes
---
# /document Workflow

You are updating documentation after code changes.

## 1. Identify Changes
- Check git diff or recent commits for modified files
- Identify which features/modules were changed
- Note any new files, deleted files, or renamed files

## 2. Verify Current Implementation
**CRITICAL**: DO NOT trust existing documentation. Read the actual code.

For each changed file:
- Read the current implementation
- Understand actual behavior (not documented behavior)
- Note any discrepancies with existing docs

## 3. Update Relevant Documentation

- **README.md**: Ensure architectural, structural, or high-level project workflow changes are accurately reflected here.
- **CHANGELOG.md**: 
  - **CRITICAL**: Before modifying the CHANGELOG, ask the user: *"Are we bumping the version? If so, what is the new semantic version and today's date?"*
  - Wait for the user's response.
  - **If the user provides a version and date:** Rename the `## [Unreleased]` header to `## [Version] - Date` and create a fresh `## [Unreleased]` block above it.
  - **If the user says no/unreleased:** Add entries under the existing `## [Unreleased]` section.
  - Use categories: Added, Changed, Fixed, Security, Removed.
  - Be concise, user-facing language.

## 4. Documentation Style Rules

✅ **Concise** - Sacrifice grammar for brevity
✅ **Practical** - Examples over theory
✅ **Accurate** - Code verified, not assumed
✅ **Current** - Matches actual implementation

❌ No enterprise fluff
❌ No outdated information
❌ No assumptions without verification

## 5. Ask if Uncertain

If you're unsure about intent behind a change or user-facing impact, **ask the user** - don't guess.