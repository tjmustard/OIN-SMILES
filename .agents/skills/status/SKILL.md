---
name: status
description: Immediately outputs the current, fully updated Living Master Plan including project objectives, architecture hypotheses, actionable steps, constraints, and a graveyard of failed approaches. Use at any point to get a snapshot of the project state.
---

# Status

This skill immediately outputs the current Living Master Plan — the complete, up-to-date snapshot of project state maintained by the Expert Co-Researcher persona.

## When to use this skill

- When the user wants a quick summary of where the project stands.
- When the user explicitly runs `/status`.
- At the start of a new session to re-orient before diving in.
- When context has been lost and the user needs to re-establish shared understanding.

## How to use it

Read the following files to compile the current state:
- `.agents/memory/activeContext.md` — current focus and next steps
- `.agents/memory/productContext.md` — product goals and decisions
- `.agents/memory/systemPatterns.md` — architectural patterns
- `spec/compiled/architecture.yml` — hypergraph state
- `spec/compiled/SuperPRD.md` — system of record

Then output the Living Master Plan immediately, covering all sections without compressing or summarizing older milestones:

### Project Objectives
The fundamental goals of the current project — what success looks like.

### Current Hypotheses & Architecture
The theoretical foundation and structural plans currently in play. What architectural bets have been made and why.

### Actionable Steps
A high-level roadmap where any immediate tasks are broken down into clear atomic steps. Flag what is blocked, in-progress, or ready to start.

### Known Constraints & Open Questions
Limitations, barriers, and uncertainties that affect execution. Include both technical constraints and product unknowns.

### The Graveyard
A strict, exhaustively detailed log of:
- Failed code iterations and why they failed
- Discarded architectural approaches and what ruled them out
- Disproven hypotheses and the evidence against them

Do not compress this section. The Graveyard prevents repeating known failures.
