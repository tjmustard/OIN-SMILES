# Example Commit Messages

This file shows real commit messages from the Hypergraph Coding Agent Framework repo, demonstrating the style conventions Claude adapts for auto-suggesting commit messages.

---

## Pattern 1: Release (Version Bump)

```
Release v0.4.1: Framework upgrades with /hyper-update and /hyper-recover
```

**When to use:** When promoting `## [Unreleased]` to a versioned release in CHANGELOG.md.  
**Format:** `Release vX.Y.Z: <brief phrase summarizing batch>`  
**Style notes:**
- Short, clear verb in present tense
- Version number prominent
- No period at end

---

## Pattern 2: Feature Addition (with multi-line body)

```
feat: add hyper-update skill to architecture and commands

Introduces /hyper-update as a smart framework upgrade mechanism. Users can
now safely upgrade the Hypergraph Coding Agent Framework while preserving
customizations to CLAUDE.md, AGENTS.md, IDE rules, and other sensitive files.

New hypergraph node: hyper_update (Module dimension)
- Auto-updates non-sensitive framework files
- Guides section-by-section collaborative merge for customizable files
- Creates timestamped backups before any mutation
- Implements approval gates to prevent data loss

Sensitive files subject to smart merge: CLAUDE.md, AGENTS.md, GEMINI.md,
IDE command bridges, and all rules files across .claude/, .clinerules/,
.roo/, .windsurf/, .cursor/ directories.

Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>
```

**When to use:** When adding a new skill or significant feature.  
**Format:** `feat: <concise subject> [+ multi-paragraph body]`  
**Style notes:**
- Conventional Commits prefix (`feat:`)
- Subject ≤50 chars
- Blank line between subject and body
- Body explains WHY, not what
- Co-Authored-By footer optional

---

## Pattern 3: Simple Addition (single-line imperative)

```
Add hyper-publish skill for timestamp-safe git publishing
```

**When to use:** For straightforward additions without complex context.  
**Format:** `Add <what> for/to <reason>`  
**Style notes:**
- Bare imperative (no `feat:` prefix)
- ≤72 chars
- Clear and direct

---

## Pattern 4: Update/Documentation

```
Update docs/ to explicitly reference /hyper-execute and fix stale paths
```

**When to use:** When updating documentation or fixing references.  
**Format:** `Update <file/area> to <reason>`

---

## Pattern 5: Fix/Correction

```
Apply linter formatting to hyper-document SKILL.md
```

**When to use:** For corrections, formatting, or applying fixes.  
**Format:** `Apply <action> to <target>`

---

## Auto-Suggestion Rules

When `/hyper-publish` is invoked with no commit message argument, the `generate_commit_message.py` script:

1. Reads the `## [Unreleased]` block in CHANGELOG.md
2. Infers the dominant change type (`Added`, `Changed`, `Removed`)
3. Examines recent commits (last 5) for style hints
4. Generates a subject line ≤72 chars using one of the patterns above
5. Optionally adds a body if multiple bullets are present in CHANGELOG

**Fallback:** If CHANGELOG `[Unreleased]` is empty, suggests a generic message based on changed file paths.
