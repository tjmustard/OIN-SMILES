---
description: Walk the user through a practical example of the Hypergraph Coding Agent Framework
---

# /tutorial Workflow

This workflow guides the user through the **Hypergraph Coding Agent Framework: A Practical Walkthrough** using the example from `docs/Tutorial.md`.

The tutorial demonstrates how to add an "Email Newsletter Subscription" feature to an existing system, illustrating the conversational specification engine and the execution engine.

## Agent Instructions

When the user runs the `/tutorial` command, follow these steps to walk them through the process. Present the information sequentially, waiting for the user to acknowledge or participate in each step before moving on. Keep your tone helpful and instructive.

### Step 1: Introduction to the Scenario
Introduce the scenario: The user has an existing project with a basic frontend and user database. Their goal is to add a form where users can subscribe to a newsletter.

Explain that the tutorial is divided into two phases: **Phase 1: The Specification Engine** and **Phase 2: The Execution Engine**.

Ask the user if they are ready to begin Step 1 of Phase 1. Wait for their response.

### Step 2: Phase 1 - Step 1 (/architect)
Explain the first step: Using the `/architect` command.
In a real scenario, the user would run `/architect I want to add an email newsletter subscription box to the homepage.`

Demonstrate how the Architect agent would respond with a paced interview:
- *Agent Question 1:* What is the primary data output? Are we storing these emails in our own database, or pushing them directly to a third-party API like Mailchimp?
- *Agent Question 2:* Are there any restrictions on who can subscribe?
- *User Response:* 1. Store in our own Postgres database. 2. Open to the public.

Explain that this interview continues until the agent creates `spec/active/Draft_PRD.md`.
Ask the user if they are ready for Step 2. Wait for their response.

### Step 3: Phase 1 - Step 2 (/redteam)
Explain the second step: Using the `/redteam` command.
**Emphasize:** The user must open a **fresh chat context** so the Red Team doesn't get confused by the Architect's conversation.

The user runs `/redteam Analyze the Draft_PRD.md`.
The Red Team reads the draft and `architecture.yml` and outputs `spec/active/RedTeam_Report.md`.

Show a sample Red Team excerpt:
- *Mutation Conflict:* No rate limiting on public submission may cause a DoS attack.
- *Missing NFRs:* No Double Opt-in (verification email) specified.

Ask the user if they are ready for Step 3. Wait for their response.

### Step 4: Phase 1 - Step 3 (/resolve)
Explain the third step: Using the `/resolve` command.
**Emphasize:** Open another **fresh chat context**.

The user runs `/resolve`. The Resolution agent presents forced trade-offs based on the Red Team Report.

Show a sample Resolution interaction:
- *Agent:* Red Team found bot submission vulnerabilities.
- *Option A:* Google reCAPTCHA v3 & IP-based rate limit (High effort, maximum security).
- *Option B:* Hidden honeypot field & DB uniqueness constraint (Low effort, moderate security).
- *User Response:* Option B for now.

Explain that the agent then compiles the final `SuperPRD` and `MiniPRD_Newsletter.md` into `spec/compiled/`. The user must then run the archival script to clear drafts: `python .agents/scripts/archive_specs.py Newsletter_Feature`.

Ask the user if they are ready to move to Phase 2. Wait for their response.

### Step 5: Phase 2 - Step 4 (The Builder)
Explain the Execution Engine. **Emphasize:** Open a **fresh chat context**.

The user asks their standard coding agent to implement the chosen MiniPRD.
*Prompt:* "Implement spec/compiled/MiniPRD_Newsletter.md. When finished, execute python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [node_ids_you_changed]"

Explain that the agent writes the code, runs the script, and the script flags nodes in `architecture.yml` as `status: needs_review`.

Ask the user if they are ready for the final step. Wait for their response.

### Step 6: Phase 2 - Step 5 (/audit)
Explain the final step: Using the `/audit` command.
**Emphasize:** Wait for the Builder to finish, then open a **fresh chat context**.

The user runs `/audit spec/compiled/MiniPRD_Newsletter.md`.
The Auditor reads the code, verifies implementation (like the honeypot), and checks the `architecture.yml`. It rewrites the YAML definitions to include the new POST data flow and sets the status back to `clean`.

Conclude the tutorial by summarizing that the user has successfully built a secure, audited feature without writing a single line of code themselves. Ask if they have any questions about the tutorial scenario.
