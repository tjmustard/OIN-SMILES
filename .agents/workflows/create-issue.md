---
description: Create a GitHub issue for a bug, feature, or improvement
---
# /create-issue Workflow

User is mid-development and thought of a bug/feature/improvement. Capture it fast so they can keep working.

## Your Goal

Create a complete issue with:
- Clear title
- TL;DR of what this is about
- Current state vs expected outcome
- Relevant files that need touching
- Risk/notes if applicable
- Proper type/priority/effort labels

## How to Get There

**Ask questions** to fill gaps - be concise, respect the user's time. They're mid-flow and want to capture this quickly. Usually need:
- What's the issue/feature
- Current behavior vs desired behavior
- Type (bug/feature/improvement) and priority if not obvious

Keep questions brief. One message with 2-3 targeted questions beats multiple back-and-forths.

**Search for context** systematically:
- ALWAYS check `spec/compiled/architecture.yml` (the hypergraph) to identify the specific `node_id`s that this issue will affect. 
- Use `grep_search` to find the relevant implementation files linked to those nodes.
- Note any risks, dependencies, or architectural constraints you spot.

**Do not guess** - Use First Principles. If the system architecture impact isn't obvious, state that in the issue.

**Keep it fast** - Total exchange under 2min. Be conversational but brief. Get what you need, create ticket, done.

## Behavior Rules

- Be conversational - ask what makes sense, not a checklist
- Default priority: normal, effort: medium (ask only if unclear)
- Max 3 files in context - most relevant only
- Bullet points over paragraphs