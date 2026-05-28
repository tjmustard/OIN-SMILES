---
name: deepdive
description: Temporarily pauses the main development plan to exhaustively explore, research, or brainstorm a specific topic using First Principles thinking and strict epistemic humility. Use when the user needs to deeply understand a concept before making a decision.
---

# Deep Dive

This skill temporarily pauses the broader development plan and exhaustively researches or brainstorms a specific topic. It applies First Principles thinking and enforces strict epistemic humility — no hallucinated data or literature.

## When to use this skill

- When the user wants to deeply understand a specific concept, technology, or design decision before proceeding.
- When the user specifies `/hyper-deepdive [topic]`.
- When an architectural decision requires thorough research before committing.

## How to use it

1. **Pause the Living Master Plan**
   - Acknowledge that you are temporarily pausing plan updates.
   - Note the current state so you can return to it after the dive.

2. **Exhaustively Explore the Topic**
   - Deconstruct the topic using First Principles thinking.
   - Research from multiple angles: what it is, how it works, its trade-offs, failure modes, and alternatives.

3. **Enforce Strict Epistemic Humility**
   - State your confidence level explicitly for specific mechanisms, papers, or technical constraints.
   - If you reach the edge of your verifiable knowledge, state this clearly and ask the user to provide data or context.
   - **Do not hallucinate literature or data.**
   - Provide standard reference formats (e.g., DOIs, PubMed IDs) if referencing scientific literature.

4. **Structure the Output**
   Present findings clearly with:
   - **Core Concept**: What it is, the problem it solves
   - **Mechanics**: How it works underneath
   - **Trade-offs**: Why you'd choose this vs. alternatives
   - **Failure Modes**: What breaks and how to debug
   - **Relevance**: How this applies to the current project/decision

5. **Integrate Findings**
   - Once the deep dive is complete, ask the user: "Would you like to integrate any of these findings into the Living Master Plan or the current specification?"
   - If yes, update the relevant memory or spec files accordingly.
