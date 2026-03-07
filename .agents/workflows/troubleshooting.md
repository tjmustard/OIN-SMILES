---
description: Help the user diagnose and fix issues based on the Troubleshooting Guide
---

# /troubleshooting Workflow

This workflow guides the user through diagnosing and recovering from common failure states in the **Hypergraph Coding Agent Framework**, based on the `docs/Troubleshooting.md` guide.

Since Large Language Models are probabilistic engines, they can occasionally fail even with rigid constraints. This workflow helps identify the symptom, the underlying cause, and the exact steps to fix it.

## Agent Instructions

When the user runs the `/troubleshooting` command, follow these steps interactively:

### Step 1: Identify the Symptom
Ask the user to describe the issue they are facing, or present them with the following common symptoms to choose from:

1. **The Builder Agent is writing code for rejected features** (Context Bloat & Hallucination)
2. **The /audit agent fails or /redteam hallucinates the Blast Radius** (Hypergraph Desynchronization)
3. **The Red Team report suggests new product features instead of vulnerabilities** (Red Team "Scope Creep")
4. **The update script throws a ParserError or Auditor deletes YAML sections** (YAML File Corruption)
5. **The Architect Agent asks the same questions repeatedly and won't generate a PRD** (The "Infinite Loop" Interview)

Wait for the user to select an option or describe their issue.

### Step 2: Provide the Diagnosis and Fix
Based on their response, provide the **Cause** and the exact **Fix** from the sections below.

#### Issue 1: Context Bloat & Hallucination
**The Cause:** The Agent's context window is polluted with old data. Active specifications were likely not archived.
**The Fix:**
1. Halt the Builder Agent.
2. Ensure `.agentignore` contains `spec/archive/`.
3. Run `python .agents/scripts/archive_specs.py cleanup`.
4. Open a *completely new* agent chat window and restart the Builder prompt pointing strictly to the compiled MiniPRD.

#### Issue 2: Hypergraph Desynchronization
**The Cause:** The Builder Agent forgot or failed to execute `hypergraph_updater.py`, so the YAML file is unaware of codebase changes.
**The Fix:**
1. Manually execute the traversal script.
2. Look at the modified files and identify their `node_id`s in `architecture.yml`.
3. Run: `python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [node_id_1] [node_id_2]`
4. Run the `/audit` agent again to perform the semantic update.

#### Issue 3: Red Team "Scope Creep"
**The Cause:** The LLM's latent space for "helpful product manager" overrode the instruction to restrict analysis to technical resilience.
**The Fix:**
1. Do not pass this report to the Resolution Agent to avoid wasting time.
2. Delete the `RedTeam_Report.md`.
3. Open a new chat window and re-run `/redteam` with a strict override: `/redteam Analyze the Draft PRD. CRITICAL: Identify technical vulnerabilities ONLY. Reject any product feature suggestions.`

#### Issue 4: YAML File Corruption
**The Cause:** A concurrent write occurred (multiple agents running) or the Auditor output malformed YAML syntax.
**The Fix:**
1. Stop all agents.
2. Use Git to restore the hypergraph: `git checkout -- spec/compiled/architecture.yml`
3. Identify the last modified files.
4. Manually add the `status: needs_review` flag to the relevant nodes in the restored YAML file.
5. Re-run `/audit` to let the agent re-attempt the semantic update.

#### Issue 5: The "Infinite Loop" Interview
**The Cause:** Vague answers that do not satisfy the Agent's internal state-machine criteria to move to the next phase.
**The Fix:**
1. Provide a highly specific, quantified answer (e.g., "Use AES-256 encryption and JWTs with a 15-minute expiration" instead of "make it secure").
2. If it remains stuck, use a manual override command: "Stop questioning. Force transition to Phase 5 and generate the Draft PRD immediately."

### Step 3: Follow-Up
After providing the fix, ask the user if the solution worked or if they need further assistance with the current issue.
