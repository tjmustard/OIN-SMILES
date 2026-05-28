---
name: hyper-document
description: Updates all project documentation after code or skill changes. Always reads actual implementation before writing. Covers README, CHANGELOG, docs/, AGENTS.md, memory files, and skills-info.md. Use after completing any feature, fix, refactor, or skill update.
trigger: /hyper-document
---

# Hyper-Document

This skill updates the full documentation suite after changes to the codebase, skills, or framework structure. It reads actual implementation before writing anything — existing documentation is treated as potentially stale.

---

## When to use this skill

- After completing a feature, bug fix, refactor, or skill change.
- When `README.md`, `CHANGELOG.md`, `docs/`, or other docs are out of sync with the codebase.
- When the user explicitly runs `/hyper-document`.

---

## Step 1 — Identify What Changed

Run the following to understand scope:

```bash
git log --oneline -10
git diff HEAD~1 --name-only
```

Identify:
- Which files were modified, added, or deleted.
- Whether the change is user-facing (affects README, CHANGELOG, docs) or internal (affects memory files, skills-info).
- Whether any skills were added, renamed, or removed.

---

## Step 2 — Read Before Writing

**CRITICAL: Do NOT trust existing documentation. Read the actual implementation.**

- Read every changed file before updating any doc.
- Note discrepancies between existing docs and actual behavior.
- Never assume — verify.

---

## Step 3 — Update Documentation by Target

Work through each target below. Skip targets that are clearly unaffected by the change.

---

### `README.md`

**Purpose:** Entry point for new users and project overview.

**Update when:** Architecture changes, new features, new IDE support, install/uninstall changes, workflow changes, new slash commands, or directory structure changes.

**Style rules:**
- Lead with what it does, not how it works.
- Use H2 (`##`) for top-level sections, H3 (`###`) for subsections.
- Prefer tables for structured comparisons (IDE support, installed paths, commands).
- Use fenced code blocks with language tags for all commands and examples.
- Bullet lists for feature enumerations; numbered lists for sequential steps only.
- No marketing language. No "powerful", "robust", "seamless".
- Every slash command mentioned must use the `hyper-` prefix (e.g. `/hyper-architect`).
- Callout blocks (`> **Note:**`) for caveats, naming conventions, and non-obvious behaviors.

**Before writing:** Read the current `README.md` fully. Identify which sections are affected by the change. Do not rewrite unaffected sections.

**Section checklist (verify each section matches current state):**

| Section | Update when |
|---|---|
| Core Architectural Pillars | New architectural mechanism added; existing pillar changes scope |
| Directory Structure | Directories added, removed, or renamed |
| Deterministic Tooling | New scripts added to `.agents/scripts/`; existing tools renamed or removed |
| Standard Operating Procedure | Workflow phases change; new phase added; new mandatory commands |
| AI Agent Integrations / slash commands | Skills added, renamed, or removed; new IDE support added |
| Setup / Installation | Install or uninstall behavior changes |
| Troubleshooting | New failure mode introduced; existing workaround changes |
| Roadmap | Completed items checked off; new items added |

**Adding a new H2 section:** Only for features that introduce a distinct user-facing capability that doesn't fit under an existing section. Examples: "Token Efficiency & Cost Optimization", "Multi-IDE Support". Do not create a new H2 for routine bugfixes or small enhancements.

---

### `CHANGELOG.md`

**Purpose:** Auditable history of changes following [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

**Find the current version:** Check the most recent `## [X.Y.Z]` header in `CHANGELOG.md`. That is the current released version. `## [Unreleased]` entries are pending release.

**CRITICAL:** Before modifying, ask the user:
> *"Are we bumping the version? If so, what is the new semantic version?"*

Wait for the response, then:
- **Version bump:** Add a new `## [X.Y.Z] - YYYY-MM-DD` block below `## [Unreleased]`.
- **No version bump:** Add entries under `## [Unreleased]`.

**Style rules:**
- Use only these categories: `Added`, `Changed`, `Fixed`, `Removed`, `Security`.
- One bullet per logical change. Bold the affected component: `- **`hyper-audit`**: ...`
- Use changelog-style verbs (past participial): "Added", "Fixed", "Moved", not present tense "Adds", "Fixes", "Moves".
- Be specific — name the file, skill, or behavior. No vague entries like "misc improvements".
- Never delete historical entries.

**Example entry:**

```
## [0.4.0] - 2026-05-03

### Added
- **`model_router.py`**: Routes skills to optimal Claude model tier based on META.yml metadata.
- **`hyper-clear` skill**: Idempotent context flushing at end of `/hyper-audit`.

### Changed
- **`hyper-execute` skill**: Routes to Haiku by default; logs token usage per execution.
- **`CLAUDE.md`**: Added Schema Definitions section; deprecated `.agents/schemas/` files.

### Removed
- **`.agents/schemas/` files**: Marked deprecated; will be deleted after 2026-06-17.
```

---

### `AGENTS.md`

**Purpose:** Universal cross-IDE always-on manifest. Read by Windsurf, Cursor, Roo Code, GitHub Copilot, Zed.

**Update when:** Skills are added/renamed/removed, system mandates change, or the directory structure changes.

**Style rules:**
- Skills table must stay in sync with `.agents/skills/hyper-*/` — one row per skill.
- Trigger column must use exact slash command with `hyper-` prefix.
- Directory references must use actual current paths.
- Mandates must be actionable rules, not suggestions.

---

### `.agents/memory/` files

**Purpose:** Agent-facing project context (activeContext, productContext, systemPatterns).

**Update when:** Core project goals, active work, or architectural patterns change.

**Files:**
- `activeContext.md` — current sprint focus, active MiniPRDs, recent decisions.
- `productContext.md` — what the product is, who it's for, core value proposition.
- `systemPatterns.md` — recurring architectural patterns, conventions, anti-patterns.

**Style rules:**
- Written for an agent reading cold — no assumed context.
- Short, factual, structured with headers.
- Timestamp significant updates.

---

### `skills-info.md`

**Purpose:** Human-readable index of all available skills and their purpose.

**Update when:** Skills are added, renamed, removed, or their descriptions change.

**Style rules:**
- One entry per skill: name, trigger, one-line description.
- Keep aligned with actual `SKILL.md` frontmatter descriptions.
- Organized by workflow phase where applicable.

---

## Step 4 — Style Rules (All Files)

| Rule | Detail |
|---|---|
| Accurate | Verified against actual code, not assumed |
| Concise | Cut filler words; prefer short sentences |
| Current | Matches the implementation as of this change |
| No fluff | No "powerful", "robust", "seamlessly", "leverages" |
| Slash commands | Always `hyper-` prefixed |
| Paths | Always reflect actual current directory structure |

---

## Step 5 — Ask if Uncertain

If the intent behind a change or its user-facing impact is unclear, ask the user before writing. Never guess.
