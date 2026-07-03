---
name: hyper-tutorial-audit
description: Run this after /hyper-tutorial-run to turn tutorial failures into a structured, constraint-aware goal prompt. Interviews the user about their goal and what may/may not be changed, then generates a paste-ready prompt for a new Claude Code fix session. Does NOT re-run tutorial commands.
trigger: /hyper-tutorial-audit
---

# Tutorial Audit

This skill is designed to run **after** `/hyper-tutorial-run` in the same conversation. It reads the tutorial, reviews what failed during the walk-through, interviews the user about their goal and constraints, and generates a structured goal prompt that can be pasted at the start of a fresh Claude Code session to drive an automated fix loop.

> **Do not re-run tutorial commands.** This skill is read-only and interview-driven. The fix work happens in the new session the goal prompt creates.

---

## Step 1 — Locate Tutorial

- If an argument was passed (tutorial name or path), resolve `tutorials/<name>/tutorial.md`.
- Otherwise, scan the current conversation for mentions of a tutorial name or path from a recent `/hyper-tutorial-run` session and use that.
- If still ambiguous, list all `tutorials/*/tutorial.md` files found in the project. Use **AskUserQuestion** to let the user select one.

---

## Step 2 — Read Tutorial

Read the full `tutorial.md`. Extract and hold in memory:
- Title (from the `#` heading)
- Overview / learning objective (the first paragraph or `## Overview` section)
- All section headings (`##`) with their ordinal numbers
- Any prerequisite declarations (package names, versions, tools)

---

## Step 3 — Review Conversation Context

Scan the current conversation (everything visible, not just recent messages) for evidence from the `/hyper-tutorial-run` session:
- Which sections the user reached or attempted
- Error messages, stack traces, or failure descriptions that appeared in chat
- Steps the user skipped, flagged as broken, or noted needed help

Build an internal summary of failures — **do not ask the user to re-describe or re-paste anything already visible in context.**

---

## Step 4 — Interview: Goal (HITL Gate #1)

Use **AskUserQuestion** with at most 2 questions:

**Question 1 — What are you trying to accomplish?** (multi-select)
- Learn the concept the tutorial teaches
- Get the feature working in my project
- Use as reference for a future task
- Reproduce a bug or a demo

**Question 2 — Which sections need attention?**

Show a checklist of all section headings from the tutorial. Pre-select any sections that appeared to have failures in Step 3, and ask the user to confirm or adjust the selection.

---

## Step 5 — Interview: Constraints (HITL Gate #2)

Use **AskUserQuestion** with at most 2 questions:

**Question 1 — What is the fix agent allowed to modify?** (multi-select)
- Tutorial text (fix outdated instructions, wrong version numbers, renamed APIs)
- Source code in the project
- Dependency versions (pyproject.toml / requirements.txt)
- Configuration files (env vars, paths, settings files)

**Question 2 — Are there things that must NOT be changed?** (single-select)
- Don't touch the source code
- Don't change the tutorial structure or learning objectives
- Don't upgrade pinned dependencies
- No hard constraints — fix whatever is needed

---

## Step 6 — Interview: Environment (HITL Gate #3 — conditional)

Only ask this gate if no environment details are visible in the conversation from Step 3. If the environment and errors are already captured, skip this gate entirely.

Use **AskUserQuestion**:

**Question 1 — What Python / tool version are you running?** (freeform)

**Question 2 — Any error messages to include that are not already in the chat?** (freeform, optional)

---

## Step 7 — Generate Goal Prompt

Using everything gathered, compose the goal prompt. Write it as a fenced markdown block so the user can copy the full prompt verbatim.

Use this structure — fill every section with real content, not placeholder text:

```
## Goal: Make "[Tutorial Title]" pass end-to-end

### Tutorial
Path: tutorials/<name>/tutorial.md

### What I'm trying to accomplish
<User's stated goals from Step 4, written as 1–2 sentences>

### Tutorial overview
<1–2 sentence learning objective extracted from the tutorial's Overview section>

### Sections that need attention
<Numbered list of sections the user flagged or that had visible failures>

### Known failures
<Failures found in conversation context — for each: section name, step or command description, error summary>

### Environment
<Python/tool version and any relevant installed packages, from conversation context or Step 6>

### What the fix agent MAY change
<Bullet list of allowed modifications from Step 5>

### What the fix agent must NOT change
<Hard constraints from Step 5, or "None" if unrestricted>

### Instructions for the fix agent

Work through `tutorials/<name>/tutorial.md`, focusing on the sections listed above:

1. **Fix in order** — address each failing section before moving to the next
2. **Dependency drift** — if a package version is wrong or missing AND dependency changes are allowed: update `pyproject.toml` (`uv add <pkg>`) AND update the corresponding tutorial text to show the correct version
3. **Tutorial text fixes** — if an instruction is outdated (wrong path, renamed API, removed option) AND tutorial text changes are allowed: update the tutorial text to match current reality
4. **Source code fixes** — only if source code changes are allowed per the constraints above
5. **Verify each fix** — after each change, re-run the failing step to confirm it passes before moving to the next
6. **Final check** — walk all sections from Section 1 to confirm no regressions introduced

### Constraints
- Do not change the tutorial's learning objectives, structure, or section order unless explicitly allowed above
- Prefer the smallest change that makes the step pass
- Keep all existing project tests passing (`uv run pytest` must exit 0) if the project has a test suite
- Never force-install packages that conflict with existing dependencies
```

---

## Step 8 — Present Output and Save Offer

Display the generated goal prompt in a fenced code block so the user can copy the full text.

Then print below the code block:

> **Next step:** Start a new Claude Code session and paste the goal prompt above as your opening message. The agent will work through the tutorial failures iteratively, respecting the constraints you've set.

Then use **AskUserQuestion** to offer a save option:
- **Option A** — Save to `tutorials/<name>/audit_<YYYY-MM-DD>.md` (use today's date)
- **Option B** — Display only — I'll copy it myself

If Option A is chosen, write the goal prompt (without the surrounding fenced code block) to that path and confirm with the filename.

---

## Constraints

- **Never run tutorial commands** — this skill is read-only and produces a goal prompt only.
- **Never modify tutorial files, source code, or any project file** during the audit (Steps 1–7 are read + interview only). Only Step 8 Option A writes a file.
- **Max 2 questions per AskUserQuestion call** — keep the interview conversational.
- **Pre-fill from context** — if conversation context from Step 3 already answers a question, pre-select or skip that question rather than asking twice.
- **Use multi-select** for all AskUserQuestion calls where choices are not mutually exclusive.
