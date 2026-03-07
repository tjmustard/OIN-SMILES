---
name: tutorial
description: Walks the user through a practical, step-by-step example of the Hypergraph Coding Agent Framework using the Email Newsletter Subscription scenario. Use when onboarding new users or demonstrating the framework workflow.
---

# Tutorial

This skill walks the user through a practical demonstration of the Hypergraph Coding Agent Framework using a real-world scenario: adding an Email Newsletter Subscription feature to an existing system.

## When to use this skill

- When onboarding a new user to the Hypergraph framework.
- When the user explicitly runs `/tutorial`.
- When the user wants a concrete walkthrough before starting their own project.

## How to use it

Present each step sequentially. Wait for the user to acknowledge or participate before moving to the next step.

### Step 1: Introduce the Scenario
Introduce the tutorial: The user has an existing project with a basic frontend and user database. Their goal is to add a newsletter subscription form.

Explain the two phases: **Phase 1: Specification Engine** and **Phase 2: Execution Engine**.

Ask: "Ready to begin Phase 1?" Wait for their response.

### Step 2: Phase 1 — `/architect`
Demonstrate the Architect's paced interview:
- *Agent:* "What is the primary data output — own database or third-party API like Mailchimp?"
- *Agent:* "Are there restrictions on who can subscribe?"
- *User:* "1. Postgres database. 2. Open to the public."

Explain the interview continues until `spec/active/Draft_PRD.md` is generated.
Ask: "Ready for Step 2?" Wait.

### Step 3: Phase 1 — `/redteam`
**Emphasize:** The user must open a **fresh chat context** to prevent contamination.

The Red Team reads the draft and `architecture.yml`, then outputs `spec/active/RedTeam_Report.md`.

Show a sample excerpt:
- *Mutation Conflict:* No rate limiting on the public submission endpoint — risk of DoS attack.
- *Missing NFR:* No Double Opt-in (verification email) — violates anti-spam compliance.

Ask: "Ready for Step 3?" Wait.

### Step 4: Phase 1 — `/resolve`
**Emphasize:** Open another **fresh chat context**.

Demonstrate forced trade-offs:
- *Option A:* Google reCAPTCHA v3 + IP-based rate limit (High effort, maximum security).
- *Option B:* Hidden honeypot field + DB uniqueness constraint (Low effort, moderate security).
- *User:* "Option B for now."

The agent compiles `SuperPRD.md` and `MiniPRD_Newsletter.md` into `spec/compiled/`. Then automatically runs: `python .agents/scripts/archive_specs.py Newsletter_Feature`.

Ask: "Ready for Phase 2?" Wait.

### Step 5: Phase 2 — The Builder
**Emphasize:** Open a **fresh chat context**.

The Builder is prompted: "Implement `spec/compiled/MiniPRD_Newsletter.md`. When finished, run `python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [node_ids]`."

The script flags the modified nodes in `architecture.yml` as `status: needs_review`.

Ask: "Ready for the final step?" Wait.

### Step 6: Phase 2 — `/audit`
**Emphasize:** Wait for the Builder to fully finish, then open a **fresh chat context**.

The Auditor reads the code, verifies the honeypot was implemented correctly, and checks the deterministic tests. It then rewrites the `architecture.yml` definitions and sets node status back to `clean`.

Conclude: "Feature complete. You successfully built a secure, audited feature without writing a single line of code yourself." Ask if they have questions.
