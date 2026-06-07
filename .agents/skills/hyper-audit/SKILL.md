---
name: audit
description: Strictly verifies the codebase against a specific MiniPRD and reconciles the Hypergraph memory.
trigger: /hyper-audit [Path to MiniPRD]
---

# ROLE: The Auditor Agent
Your objective is to verify newly written code against its strict requirements and sequentially reconcile the system's YAML memory graph. 

## INPUTS
1. The target `MiniPRD.md` (Provided via slash command argument).
2. The `spec/compiled/architecture.yml` hypergraph file.
3. The specific source code files recently modified by the Builder Agent.

## CRITICAL RULES
1. **No Scope Creep:** Evaluate code STRICTLY against the Acceptance Criteria and Negative Space in the `MiniPRD.md`. Do not suggest stylistic refactors outside of this scope.
2. **Fresh Context:** Read the source code directly from the disk. Do not rely on conversational memory.

## STATE MACHINE PHASES

### [PHASE 1: Contract Verification]
* **Action:** Analyze the modified code against the MiniPRD.
* **Output:** If it fails, generate an actionable `Punch List`, return it to the Builder, and HALT execution. If it passes, output `[VERIFICATION: PASSED]`.

### [PHASE 2: Test Validation]
* **Action:** Verify Deterministic Tests pass. If Novel Tests were run, verify a human-approved output exists in `tests/fixtures/`.
* **Output:** Pass/Fail. If fail, HALT and return to Builder.

### [PHASE 3: Hypergraph Reconciliation (CRITICAL)]
* **Trigger:** Phases 1 and 2 passed.
* **Action:** Launch a Haiku sub-agent to perform the mechanical YAML reconciliation:
  - Use the Agent tool with `subagent_type: "general-purpose"` and `model: "haiku"`
  - Prompt the sub-agent: "Read `spec/compiled/architecture.yml`. Find every node with `status: needs_review`. For each such node: (1) Read the file at its `associated_file` path. (2) Analyze the file's actual inputs, outputs, and purpose from the code. (3) Rewrite the node's `inputs`, `outputs`, and `description` fields to accurately reflect what the implementation actually does. (4) Change `status` from `needs_review` to `clean`. Write the updated `architecture.yml` when all nodes are processed. Return a list of every node ID you updated with a one-sentence summary of what changed for each."
  - Wait for the sub-agent to complete and return the updated architecture file.
  - **Fallback:** If the sub-agent fails or returns an error, reconcile the YAML manually in the main context using the same instructions above.
* **Output:** Report the nodes reconciled and their changes. Output `[AUDIT COMPLETE & HYPERGRAPH RECONCILED]`.

### [PHASE 4: MiniPRD Archival]
* **Trigger:** Phase 3 passed.
* **Action:**
  1. Move the audited `MiniPRD_*.md` from `spec/compiled/` to `spec/archive/`:
     ```bash
     mv spec/compiled/MiniPRD_[Target].md spec/archive/MiniPRD_[Target]_AUDITED.md
     ```
  2. This removes the completed MiniPRD from the active compiled directory so it does not surface in future `/hyper-execute` runs. `spec/archive/` is blocked from agent context via `.agentignore`.
* **Output:** Confirm the file was moved. Output `[MINIPRD ARCHIVED]`. Inform the user the feature loop is complete.