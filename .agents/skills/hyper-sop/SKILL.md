---
name: sop
description: Explains the Master SOP and guides the user to the correct framework phase. Use when a new user needs orientation or when you need to re-anchor to the correct workflow phase.
trigger: /hyper-sop
---

# /hyper-sop — Master Standard Operating Procedure Guide

This skill guides the user through the **Master Standard Operating Procedure (SOP)** of the Hypergraph Coding Agent Framework.

The framework relies on a strict sequential model — Spec-First, Deterministic Memory, and Automated Auditing — to prevent context collapse and hallucinated requirements.

## Step 1: Introduction

Briefly introduce the Hypergraph framework's core paradigm:
- The workflow is strictly sequential.
- Specialized agents handle distinct phases (Architect, Red Team, Resolution, Auditor).
- Deterministic Python scripts manage state so LLMs never traverse the graph probabilistically.
- New context windows (new conversations) prevent cross-contamination between adversarial agents.

## Step 2: Determine User State

Ask the user which phase they are in:

- **Phase -1: Legacy Onboarding** — Integrating the framework into an existing project
- **Phase 0: System Initialization** — Starting a brand new greenfield project
- **Phase 1: Specification** — Planning and designing a new feature
- **Phase 2: Execution** — Writing code based on a completed MiniPRD

## Step 3: Provide Phase-Specific Guidance

Wait for the user's response, then provide the exact steps:

**Phase -1 (Legacy Onboarding):**
1. Run `/hyper-discover` — scans the codebase and populates `spec/compiled/architecture.yml`.
2. Review the generated YAML for accuracy.
3. Run `/hyper-baseline` — generates the initial `SuperPRD.md`.
4. Once verified, proceed to Phase 1.

**Phase 0 (System Initialization):**
1. Ensure the template is cloned with `.agents/` and `spec/` directories intact.
2. Install the dependency: `pip install pyyaml`
3. Give execution permissions: `chmod +x .agents/scripts/*.py`
4. The workspace is ready. Proceed to Phase 1.

**Phase 1 (Specification):**
1. Run `/hyper-architect` — begins the requirements extraction interview (creates `Draft_PRD.md` in `spec/active/`).
2. **Start a new conversation** and run `/hyper-redteam` — adversarial analysis (creates `RedTeam_Report.md` in `spec/active/`).
3. **Start a new conversation** and run `/hyper-resolve` — resolves trade-offs and compiles `SuperPRD.md` and `MiniPRD` files into `spec/compiled/`.
4. The resolve agent will automatically run `python .agents/scripts/archive_specs.py [Feature_Name]` to flush active drafts.

**Phase 2 (Execution):**
1. **Start a new conversation** and instruct the Builder to implement a specific MiniPRD from `spec/compiled/`.
2. The Builder MUST run `python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [modified_node_ids]` after code is written.
3. **Start a new conversation** and run `/hyper-audit spec/compiled/MiniPRD_[Target].md` to verify code against the contract.

## Step 4: Remind About Testing and Iteration

- **Novel Tests:** Subjective/AI-generated outputs go to `tests/candidate_outputs/` for human review. Once approved, move to `tests/fixtures/` and update the MiniPRD test definition to `deterministic`.
- **Iterative Updates:** Future features follow the identical Phase 1 → Phase 2 loop. Agents detect the existing `architecture.yml` and perform Delta Extractions automatically.

## Step 5: Next Actions

Ask the user if they are ready to begin their chosen phase or if they have questions about any specific command.
