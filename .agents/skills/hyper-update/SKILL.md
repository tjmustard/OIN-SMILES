---
name: update
description: Interactive smart-upgrade command. Fetches latest upstream framework, auto-updates non-sensitive files, and guides user through section-by-section collaborative merge of sensitive customizable files. Backs up any replaced or merged file to .agents/.backup/YYYY-MM-DD/.
---

# ROLE: The Upgrade Agent

Your objective is to safely upgrade the local Hypergraph Coding Agent Framework installation while preserving user customizations. You act as a careful, non-destructive upgrade orchestrator.

## CRITICAL RULES

1. **Never mutate without backup.** Before replacing or merging ANY file, create `.agents/.backup/YYYY-MM-DD/` (timestamped directory) and copy the original file there first.
2. **Section-by-section approval gate.** When merging, do NOT advance to the next `###` section until the user provides explicit approval of the current section. Approval is explicit: accept upstream, keep local, or edit manually — each must be confirmed by the user.
3. **Smart file categorization.** Sensitive customizable files (CLAUDE.md, AGENTS.md, rules files) require merge/replace/skip prompts. Non-sensitive files (.agents/ skill implementations, spec/, tests/) are auto-updated silently.

---

## EXECUTION PHASES

### [PHASE 1: Preflight & Fetch]

**Objective:** Prepare the upgrade environment.

**Action:**
1. Verify git is available (`git --version`).
2. Fetch the upstream repository into a temp directory:
   ```bash
   git clone --depth=1 --branch main https://github.com/tjmustard/Hypergraph-Coding-Agent-Framework.git /tmp/hcaf-upstream-$(date +%s)
   ```
3. Report: "Upstream fetched successfully. Scanning for changes..."
4. Do NOT yet create the backup directory — only create it when the first file is actually replaced/merged.

---

### [PHASE 2: Auto-Update Non-Sensitive Files]

**Objective:** Update framework core files without user interaction.

**Action:**

Silently copy these directories from upstream, if changed:
- `.agents/skills/` (all skill SKILL.md files)
- `.agents/scripts/` (Python scripts)
- `.agents/schemas/` (schema definitions)
- `spec/` (specification templates)
- `tests/` (test fixtures, candidate outputs)
- `.agentignore`

Report: "✅ Auto-updated: .agents/skills/, .agents/scripts/, spec/, tests/, .agentignore" (no individual prompts).

---

### [PHASE 3: Smart Merge of Sensitive Files]

**Objective:** Guide user through selective upgrade of customizable files.

Sensitive files are:
- `CLAUDE.md`
- `AGENTS.md`
- `GEMINI.md`
- `.claude/commands/*.md` (all command bridge files)
- `.clinerules/*.md`
- `.roo/rules/*.md`, `.roo/rules-code/*.md`
- `.cursor/rules/*.mdc`
- `.windsurf/rules/*.md`

**For EACH sensitive file that exists BOTH locally AND in upstream:**

1. **Diff check:** Run `diff -u local-version upstream-version > /tmp/diff.txt`
   - If no diff → Report: "✅ Already up to date: FILENAME" → skip to next file
   - If diff exists → proceed to step 2

2. **Prompt:** Show the filename and a one-line summary of the change.
   ```
   📝 File: CLAUDE.md
   Changed lines: 45 (schema definitions, tool names table)
   ```
   Use **AskUserQuestion**:

   ```
   How do you want to handle [FILENAME]?

   - Option A: Replace — overwrite local with upstream version (backup created)
   - Option B: Merge — section-by-section collaborative merge
   - Option C: Skip — keep local version unchanged
   ```

3. **If Replace:**
   - Create `.agents/.backup/YYYY-MM-DD/` if not exists
   - Copy local file to `.agents/.backup/YYYY-MM-DD/FILENAME`
   - Copy upstream file to local path
   - Report: "✅ Replaced FILENAME (backup: .agents/.backup/YYYY-MM-DD/FILENAME)"

4. **If Skip:**
   - Report: "⏭️  Skipped FILENAME (keeping local version)"

5. **If Merge:**
   - Create `.agents/.backup/YYYY-MM-DD/` if not exists
   - Copy local file to `.agents/.backup/YYYY-MM-DD/FILENAME`
   - Parse BOTH local AND upstream versions into `###`-level sections
   - For EACH section that differs:
     a. Display:
        - Section heading
        - **Upstream changes:** (show the upstream version of this section)
        - **Your local version:** (show the local version of this section)
        - **Unified diff** (3 lines of context)
     b. Use **AskUserQuestion**:

        ```
        How do you want to handle this section: [SECTION_HEADING]?

        - Option A: Keep local — use your existing version of this section
        - Option B: Accept upstream — use the upstream version of this section
        - Option C: Edit manually — paste a custom merged version (you will be prompted)
        ```
     c. **CRITICAL: Do NOT advance to the next section until the user explicitly approves this one.** Wait for their response.
     d. If "Keep local" → use local version of this section
     e. If "Accept upstream" → use upstream version of this section
     f. If "Edit manually" → ask: "Paste your custom merged version below (or leave blank to skip):" and use their input
   - After all sections approved, assemble merged file and write to local path
   - Report: "✅ Merged FILENAME (backup: .agents/.backup/YYYY-MM-DD/FILENAME)"

---

### [PHASE 4: Summary & Verification]

**Objective:** Confirm all changes and clean up.

**Action:**

1. List all files processed:
   ```
   📊 Upgrade Summary
   ✅ Auto-updated: .agents/skills/, .agents/scripts/, spec/, tests/ (4 categories)
   ✅ Replaced: AGENTS.md
   ✅ Merged: CLAUDE.md (3 sections), .clinerules/hypergraph-agent.md (1 section)
   ⏭️  Skipped: GEMINI.md, .roo/rules/testing.md
   ```

2. If any files were replaced/merged, report backup location:
   ```
   🔐 Backup created: .agents/.backup/2026-05-04/
   ```

3. Report completion:
   ```
   ✅ Upgrade complete! Your customizations have been preserved.
   Tip: Commit the changes with 'git commit -m "chore: update Hypergraph framework"'
   ```

---

## Negative Space (Constraints)

- **DO NOT** auto-update IDE-specific bridge files under `.claude/commands/`, `.clinerules/`, `.roo/`, etc. without user approval — these may contain project-specific overrides.
- **DO NOT** create the backup directory until the first file is actually replaced/merged.
- **DO NOT** skip section-by-section approval in merge mode. Each `###` section requires explicit user approval before advancing.
- **DO NOT** merge files at the character level — only at `###` heading boundaries.
- **DO NOT** perform this upgrade if the repo has uncommitted changes to framework files. Warn the user first.

---

## When to use this skill

- User runs `/hyper-update` to upgrade the framework.
- The local installation is already complete (`.agents/` exists).
- The user wants to pull upstream improvements while keeping customizations.
- After running `HACF-install.sh` upgrade mode fails to merge safely.

---

## Integration & Idempotency

Running `/hyper-update` multiple times is safe:
- Files with no upstream changes are skipped.
- Files with changes trigger the merge/replace/skip flow again.
- Backups are timestamped, so multiple runs don't overwrite each other.
