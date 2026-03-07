# **Troubleshooting Guide: Agent Failure States**

Large Language Models are probabilistic engines. Even with rigid constraints, they will occasionally fail. This guide details how to diagnose and recover from the most common failure states in the Hypergraph Coding Agent Framework architecture.

## **1\. Context Bloat & Hallucination**

**The Symptom:** The Builder Agent starts writing code for features that were discussed during the Interview phase but were explicitly rejected or changed during the Resolution phase.

**The Cause:** The Agent's context window is polluted with old data. You likely forgot to archive the active specifications.

**The Fix:**

1. Halt the Builder Agent.  
2. Ensure .agentignore contains spec/archive/.  
3. Run python .agent/scripts/archive\_specs.py cleanup.  
4. Open a *completely new* agent chat window and restart the Builder prompt pointing strictly to the compiled MiniPRD.

## **2\. Hypergraph Desynchronization**

**The Symptom:** The /audit agent fails to reconcile the hypergraph, or the /redteam agent hallucinates the Blast Radius.

**The Cause:** The Builder Agent forgot or failed to execute hypergraph\_updater.py, meaning the YAML file is unaware that the codebase has changed.

**The Fix:**

1. You must manually execute the traversal script.  
2. Look at the files the Builder just modified. Identify their corresponding node\_ids in architecture.yml.  
3. Run: python .agent/scripts/hypergraph\_updater.py spec/compiled/architecture.yml \[node\_id\_1\] \[node\_id\_2\]  
4. Run the /audit agent again to perform the semantic update.

## **3\. Red Team "Scope Creep"**

**The Symptom:** The Red Team report is full of suggestions for entirely new product features (e.g., "You should add a machine learning recommendation engine to this simple to-do app") instead of focusing on technical vulnerabilities.

**The Cause:** The LLM's latent space weights for "helpful product manager" overrode the system prompt's instruction to restrict analysis to technical resilience.

**The Fix:**

1. Do not pass this report to the Resolution Agent; it will waste your time asking about trade-offs for features you don't want.  
2. Delete the RedTeam\_Report.md.  
3. Open a new chat window and re-run /redteam, but append a strict override: /redteam Analyze the Draft PRD. CRITICAL: Identify technical vulnerabilities ONLY. Reject any product feature suggestions.

## **4\. YAML File Corruption**

**The Symptom:** The hypergraph\_updater.py script throws a yaml.parser.ParserError or the Auditor agent completely deletes sections of the architecture.yml file.

**The Cause:** Either a concurrent write occurred (you ran two agents at once), or the Auditor agent output malformed YAML syntax.

**The Fix:**

1. Stop all agents.  
2. Use Git to restore the hypergraph to the last known good state:  
   git checkout \-- spec/compiled/architecture.yml  
3. Identify the last modified files.  
4. Manually add the status: needs\_review flag to the relevant nodes in the restored YAML file.  
5. Re-run /audit to let the agent re-attempt the semantic update.

## **5\. The "Infinite Loop" Interview**

**The Symptom:** The Architect Agent keeps asking the same clarifying questions over and over, never generating the Draft PRD.

**The Cause:** You are likely giving vague answers that do not satisfy the Agent's internal state-machine criteria for moving to the next phase.

**The Fix:**

1. Provide a highly specific, quantified answer. (e.g., Instead of "make it secure", say "Use AES-256 encryption and JWTs with a 15-minute expiration").  
2. If it is stuck, use a manual override command: Stop questioning. Force transition to Phase 5 and generate the Draft PRD immediately.