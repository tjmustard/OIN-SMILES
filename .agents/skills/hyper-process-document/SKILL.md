---
name: process-document
description: Documents the process, methodology, and decisions of the current session as a reproducible narrative. Saves to spec/process/ in the workspace so the journey can be understood and recreated by others.
trigger: /hyper-process-document
argument-hint: "Brief title for this process document (e.g. 'Dynamic Workflow Engine')"
---

# Hyper-Process-Document

Write a process document that captures the *how* of this session — the problem, approach, decisions, and steps taken — so that someone reading it cold can understand and reproduce the results.

This is distinct from `/hyper-handoff`, which is forward-looking (context for a new agent to continue). This document is retrospective and human-readable: institutional knowledge capture.

---

## Step 1 — Determine Output Path

Create the output directory if it does not exist:

```bash
mkdir -p spec/process
```

Name the file: `process_<YYYYMMDD_HHMMSS>_<slug>.md`

- Use the current timestamp for `<YYYYMMDD_HHMMSS>`
- Derive `<slug>` from `$ARGUMENTS` by lowercasing and replacing spaces/special characters with hyphens (e.g., "Dynamic Workflow Engine" → `dynamic-workflow-engine`)
- If no `$ARGUMENTS` provided, use `session` as the slug

---

## Step 2 — Gather Context

Before writing, read and recall the following to reconstruct an accurate picture:

```bash
git log --oneline -20
```

```bash
git diff HEAD~5..HEAD --name-only
```

- Read `.agents/memory/activeContext.md` — current sprint state and recent decisions
- List `spec/active/` — any in-progress PRDs or artifacts touched this session
- Recall the conversation: what problems were raised, what approaches were tried, what choices were made and why

---

## Step 3 — Write the Process Document

Write the document to the path determined in Step 1. Use this structure:

```markdown
# Process Document: <title from $ARGUMENTS or inferred>

**Generated:** <ISO 8601 timestamp>
**Session Focus:** <from $ARGUMENTS, or inferred from conversation context>

## Problem Statement

<What problem or goal prompted this session? What was the starting condition or gap?
1–3 sentences. Written for someone with no prior context.>

## Starting State

<Describe the codebase or system state at the beginning of the session.
Reference the git commit SHA or describe key pre-existing artifacts.
What existed, what was missing, what was broken?>

## Approach & Methodology

<2–4 sentences on the overall strategy employed.
Was this exploratory? Spec-driven? Iterative debugging?
Which HACF framework phases were invoked (e.g., /hyper-architect → /hyper-execute)?
What was the sequencing rationale?>

## Steps Taken

<Numbered chronological narrative of the work performed. For each step include:
- What was attempted
- Why (the reasoning, constraint, or observation that drove it)
- What tool, skill, or command was used
- The outcome or observation

Example format:
1. Ran `/hyper-discover` to initialize the hypergraph from the existing codebase, because no
   architecture.yml existed yet. Result: spec/compiled/architecture.yml created with 12 nodes.
2. Ran `/hyper-architect` to extract requirements — identified 3 key user stories and 2 constraints.
3. Noticed the session naming hook was outputting raw text instead of JSON; traced the issue to
   SKILL.md step 3 which lacked an envelope wrapper. Fixed by wrapping output in hookSpecificOutput.>

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| <what was decided> | <what else was considered> | <why this was selected> |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| <name> | <path> | <created / updated / deleted> |

## Results & Outcomes

<What was produced? What works now that didn't before?
Be concrete — reference file paths, test outputs, or directly observable behavior.>

## How to Reproduce

<Step-by-step instructions for someone starting from a clean slate to recreate these results.
Include:
- Prerequisite state (git branch, installed tools, environment vars if any)
- Exact skills or commands to invoke, in order
- Expected output or artifact at each step
- Any gotchas or order-dependencies>

## Patterns & Lessons

<Optional. Generalizable insights from this session:
- New coding patterns or conventions adopted
- Workflow improvements discovered
- Pitfalls to avoid next time
- Framework behaviors worth noting>
```

---

## Step 4 — Redact Sensitive Information

Before writing, scan for and remove:
- API keys, tokens, secrets (any string matching `sk-`, `ghp_`, `xoxb-`, or similar)
- Passwords or credentials
- Personally identifiable information (email addresses, phone numbers, full names of non-public individuals)

Replace redacted content with `[REDACTED]`.

---

## Step 5 — Confirm Output

After saving the file, print the full path to the user:

```
Process document saved to: spec/process/process_20260616_143022_dynamic-workflow-engine.md
```

Do NOT print the document contents to the conversation — it would pollute the context this skill is trying to preserve.
