# **Hypergraph Coding Agent Framework: A Practical Walkthrough**

This tutorial demonstrates exactly how to use the Hypergraph Coding Agent Framework template by walking through a real-world scenario: **Adding an Email Newsletter Subscription feature** to an existing system.

## **The Scenario**

Your project currently has a basic frontend and a user database. You want to add a form where users can subscribe to a newsletter.

## **Phase 1: The Specification Engine**

### **Step 1: /architect**

You open your IDE and trigger the Architect skill.

**You:** /architect I want to add an email newsletter subscription box to the homepage.

**Architect Agent:**

1. What is the primary data output? Are we storing these emails in our own database, or pushing them directly to a third-party API like Mailchimp?  
2. Are there any restrictions on who can subscribe (e.g., must be a logged-in user, or open to public visitors)?

**You:** 1\. Store them in our own Postgres database for now. 2\. Open to the public.

*(The Architect continues its paced interview, asking about spam prevention and frontend UI. Once satisfied, it generates spec/active/Draft\_PRD.md)*

### **Step 2: /redteam**

You open a **fresh chat context** so the Red Team doesn't get confused by the Architect's conversation.

**You:** /redteam Analyze the Draft\_PRD.md

*(The Red Team reads the Draft PRD and queries architecture.yml to find the Blast Radius. It outputs spec/active/RedTeam\_Report.md)*

**Red Team Output (Excerpt):**

**State Mutation Conflicts:** The Draft PRD does not specify rate limiting on the public submission endpoint. A malicious actor could execute a rapid-fire script to fill the Postgres database with garbage data, causing a denial of service.

**Missing NFRs:** There is no Double Opt-in (verification email) specified. This violates standard anti-spam compliance.

### **Step 3: /resolve**

You open another **fresh chat context**.

**You:** /resolve

**Resolution Agent:**

The Red Team identified a critical vulnerability regarding bot submissions and rate limiting.

**Option A:** Implement Google reCAPTCHA v3 on the frontend and a strict IP-based rate limit on the backend (High effort, maximum security).

**Option B:** Implement a simple hidden honeypot field on the form and a basic database uniqueness constraint (Low effort, moderate security).

Which path do we take for V1?

**You:** Option B for now.

*(The Resolution Agent compiles the final SuperPRD and MiniPRD\_Newsletter.md into spec/compiled/. It then instructs you to run python .agent/scripts/archive\_specs.py Newsletter\_Feature. You run the script in your terminal, and the active drafts are safely archived.)*

## **Phase 2: The Execution Engine**

### **Step 4: The Builder**

You open a **fresh chat context** and instantiate your standard coding agent.

**You:** Implement spec/compiled/MiniPRD\_Newsletter.md. When finished, execute python .agent/scripts/hypergraph\_updater.py spec/compiled/architecture.yml \[node\_ids\_you\_changed\]

*(The agent writes the HTML form and the backend API route. It then runs the Python script in the terminal. The script flags the database\_module and homepage\_ui nodes in architecture.yml as status: needs\_review.)*

### **Step 5: /audit**

You wait for the Builder to completely finish. You open a **fresh chat context**.

**You:** /audit spec/compiled/MiniPRD\_Newsletter.md

*(The Auditor reads the code, verifies the honeypot was implemented correctly, and checks the deterministic tests. It then reads architecture.yml, sees the needs\_review flags, and rewrites the YAML definitions to include the new POST /subscribe data flow. It sets the status back to clean.)*

**Auditor Agent:**

\[VERIFICATION: PASSED\]

\[AUDIT COMPLETE & HYPERGRAPH RECONCILED\]

**Feature Complete.** You have successfully built a secure, audited feature without writing a single line of code yourself.