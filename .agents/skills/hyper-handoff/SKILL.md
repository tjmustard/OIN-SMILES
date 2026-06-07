---
name: handoff
description: Compacts the current conversation into a handoff document so a fresh agent can continue the work. Saves to the OS temp directory, not the workspace. Pass an optional argument describing what the next session will focus on.
trigger: /hyper-handoff
argument-hint: "What will the next session be used for?"
---

# Hyper-Handoff

Write a handoff document summarising the current conversation so a fresh agent can continue the work without needing this conversation history.

---

## Step 1 — Determine Output Path

Resolve the OS temp directory by running:

```bash
python3 -c "import tempfile; print(tempfile.gettempdir())"
```

Name the file: `handoff_<YYYYMMDD_HHMMSS>.md` using the current timestamp.
Save there — **not** in the current workspace.

---

## Step 2 — Read Key Artifacts (Do Not Duplicate)

Before writing, read and note paths/URLs for any existing artifacts that capture ground truth. Reference them by path rather than reproducing their content:

- `spec/compiled/architecture.yml` — hypergraph state
- `spec/active/` — any in-progress PRDs, Red Team reports, or draft specs
- `spec/compiled/` — finalized SuperPRDs and MiniPRDs
- `.agents/memory/activeContext.md` — current sprint state and recent decisions
- Recent git log (`git log --oneline -5`) — commit history

---

## Step 3 — Write the Handoff Document

Structure the document as follows:

```markdown
# Handoff: <brief title derived from conversation context>

**Generated:** <ISO 8601 timestamp>
**Next session focus:** <from $ARGUMENTS, or "general continuation" if none>

## Context Summary

<2–4 sentences: what problem is being solved and where the work stands.
No fluff. Written for an agent reading cold.>

## What Was Decided

<Bullet list of key decisions made this session. If a decision is captured in
a PRD or commit, reference the path — do not reproduce it.>

## What Was NOT Resolved

<Bullet list of open questions, blocked items, or deferred decisions.
Be specific — vague items like "needs more work" are not useful.>

## Current State of Artifacts

| Artifact | Path | Status |
|---|---|---|
| architecture.yml | spec/compiled/architecture.yml | <clean / has dirty nodes> |
| Active specs | spec/active/ | <list files or "empty"> |
| Compiled specs | spec/compiled/ | <list key files> |
| activeContext.md | .agents/memory/activeContext.md | <one-line summary> |

## Suggested Skills for Next Session

<List the HACF skills the next agent should invoke, in recommended order,
with a one-line reason for each. Choose from:>

- `/hyper-architect` — if requirements still need extraction
- `/hyper-redteam` — if a Draft PRD exists and needs adversarial review
- `/hyper-resolve` — if a Red Team report exists and needs trade-off resolution
- `/hyper-execute` — if a MiniPRD is ready to implement
- `/hyper-audit` — if code was just written and needs verification
- `/hyper-document` — if documentation is out of sync
- `/hyper-contextualize` — if agent files may have HACF framing issues
- `/hyper-status` — to reorient before any phase
- `/hyper-troubleshooting` — if something is broken or desynchronized

## Exact Next Action

<One sentence: the single most important thing the next agent should do first.
Example: "Run /hyper-redteam on spec/active/Draft_PRD.md in a new context window.">
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
Handoff document saved to: /tmp/handoff_20260606_143022.md
```

Do NOT print the document contents to the conversation — it may be long and would pollute the context this skill is trying to compact.