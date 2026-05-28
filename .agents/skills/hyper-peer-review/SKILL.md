---
name: peer-review
description: Evaluates external peer review findings against the actual codebase, separating valid issues from misunderstandings, and produces a prioritized action plan. Use when another model or reviewer has provided feedback on the current implementation.
---

# Peer Review

This skill critically evaluates external peer review findings as the team lead — verifying each finding against the actual code before accepting or rejecting it.

## When to use this skill

- When another model, reviewer, or team member has provided feedback on the current implementation.
- When the user explicitly runs `/hyper-peer-review` and pastes review findings.
- When cross-checking AI-generated code review output before acting on it.

## How to use it

1. **Receive the Findings**
   Ask the user to paste the peer review feedback if not already provided. The reviewer has less context on this project's history and decisions than you do — evaluate accordingly.

2. **Verify Each Finding**
   For EACH finding:
   - **Check if it exists** — Read the actual code. Does this issue really exist?
   - **If it doesn't exist** — Explain clearly why (already handled, reviewer misunderstood the architecture, outdated assumption).
   - **If it does exist** — Assess severity: Critical / High / Medium / Low.

3. **Produce the Summary**
   Output three sections:

   **Valid Findings (Confirmed Issues)**
   - List each confirmed issue with its severity and the specific file/function affected.

   **Invalid Findings (Rejected with Explanation)**
   - List each rejected finding with a clear explanation of why it's incorrect or inapplicable.

   **Prioritized Action Plan**
   - Ordered list of confirmed issues to fix, from highest to lowest severity.
   - For each: the specific change needed and which file/node it affects.

## Behavior Rules

- You are the team lead — do not accept findings at face value.
- Always read the code before accepting or rejecting a finding.
- Be direct but fair. If a finding is wrong, explain why clearly.
- Cross-reference `spec/compiled/architecture.yml` for architectural context when evaluating structural findings.
