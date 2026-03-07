---
name: create-issue
description: Quickly captures a bug, feature request, or improvement as a structured GitHub issue without interrupting the user's development flow. Use when the user needs to log something fast while staying focused on current work.
---

# Create Issue

This skill captures a bug, feature request, or improvement as a complete, well-structured GitHub issue — fast. The user is mid-development and wants to capture an idea without losing their flow.

## When to use this skill

- When the user is mid-development and thinks of a bug, feature, or improvement to capture.
- When the user explicitly runs `/create-issue`.
- When the user says something like "log this as an issue" or "create a ticket for this."

## How to use it

1. **Ask Targeted Questions**
   - Be concise — the user is mid-flow. Ask only what's needed to fill gaps.
   - Usually need: what's the issue/feature, current behavior vs. desired behavior, type (bug/feature/improvement) and priority if not obvious.
   - Keep it to one message with 2–3 targeted questions max. Do not do multiple back-and-forths.

2. **Search for Context**
   - Check `spec/compiled/architecture.yml` to identify the specific `node_id`s that this issue will affect.
   - Search the codebase to find the relevant implementation files linked to those nodes.
   - Note any risks, dependencies, or architectural constraints spotted.

3. **Draft the Issue**
   Create a complete issue with:
   - **Title**: Clear, concise, action-oriented
   - **TL;DR**: One-sentence summary
   - **Current State vs. Expected Outcome**: What happens now vs. what should happen
   - **Relevant Files / Nodes**: Most relevant files and `node_id`s from the hypergraph (max 3)
   - **Risk/Notes**: Architectural constraints, dependencies, or risks (if applicable)
   - **Labels**: Type (`bug`/`feature`/`improvement`), priority (`high`/`normal`/`low`), effort (`small`/`medium`/`large`)

4. **Create the Issue**
   Use the shell to run `gh issue create` with the drafted content, or output the issue body for the user to copy if GitHub CLI is unavailable.

## Behavior Rules

- Be conversational — ask what makes sense, not a checklist.
- Default priority: `normal`, effort: `medium` (ask only if unclear).
- Max 3 files in context — most relevant only.
- Use bullet points over paragraphs.
- Do not guess architecture impact — if unclear, state that in the issue.
- Total exchange: under 2 minutes.
