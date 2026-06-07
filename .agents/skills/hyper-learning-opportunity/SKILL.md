---
name: learning-opportunity
description: Pauses development mode to teach a technical concept at three increasing levels of complexity, tailored to a technical PM audience. Use when the user wants to understand something they encountered while building.
---

# Learning Opportunity

This skill pauses development mode and shifts into teaching mode — delivering a concept at three increasing levels of complexity, peer-to-peer, without oversimplification.

## When to use this skill

- When the user encounters something they don't fully understand and wants to learn it.
- When the user explicitly runs `/hyper-learning-opportunity [concept]`.
- When a development decision reveals a knowledge gap worth addressing.

## When to use this skill

- When the user wants to understand a concept, pattern, or technology they encountered.
- When the user runs `/hyper-learning-opportunity [topic]`.

## How to use it

### Target Audience
Technical PM with mid-level engineering knowledge. Understands architecture, can read code, ships production apps. Not a senior engineer, but not a beginner. Apply the **80/20 rule** — focus on concepts that compound. Don't oversimplify, but prioritize practical understanding over academic completeness.

### Deliver Three Levels

Begin by using **AskUserQuestion** to check how deep the user wants to go:

```
What depth of explanation do you want?

- Option A: Conceptual overview (Level 1 only) — what it is and why it exists
- Option B: Practical walkthrough (Levels 1–2) — concepts + mechanics + trade-offs
- Option C: Full deep dive (all 3 levels) — including implementation details and senior-engineer perspective
```

After delivering each level, use **AskUserQuestion** before advancing:

```
Ready to continue to the next level?

- Option A: Yes, continue — move to the next level
- Option B: Explain more here — go deeper on this level before advancing
```

#### Level 1: Core Concept
- What this is and why it exists
- The problem it solves
- When you'd reach for this pattern
- How it fits into the broader architecture

#### Level 2: How It Works
- The mechanics underneath
- Key trade-offs and why this approach was chosen
- Edge cases and failure modes to watch for
- How to debug when things go wrong

#### Level 3: Deep Dive
- Implementation details that affect production behavior
- Performance implications and scaling considerations
- Related patterns and when to use alternatives
- The "senior engineer" perspective — what experienced devs know that beginners don't

### Tone
- Peer-to-peer, not teacher-to-student
- Technical but not jargon-heavy
- Use concrete examples from the current codebase where possible
- Acknowledge genuine complexity honestly: "This is tricky because..."
