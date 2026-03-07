---
description: Explain the Hypergraph Coding Agent Framework Master SOP to the user
---

# /sop Workflow

This workflow is designed to help you (the agent) guide the user through the **Master Standard Operating Procedure (SOP)** of the Hypergraph Coding Agent Framework. 

The framework relies on a strict sequential model—Spec-First, Deterministic Memory, and Automated Auditing—to prevent context collapse and hallucinated requirements.

## Agent Instructions

When the user runs the `/sop` command, follow these steps sequentially:

### Step 1: Introduction
Briefly introduce the Hypergraph framework's core paradigm to the user. Explain that the workflow is strictly sequential and relies on specialized phases, deterministic graph traversal, and separate context windows to prevent system corruption.

### Step 2: Determine User State
Ask the user to identify their current development phase so you can provide relevant instructions. Give them the following options:
- **Phase -1: Legacy Onboarding** (Integrating the framework into an existing project)
- **Phase 0: System Initialization** (Starting a brand new greenfield project)
- **Phase 1: Specification** (Planning and designing a new feature)
- **Phase 2: Execution** (Writing code based on a completed MiniPRD)

### Step 3: Provide Phase-Specific Guidance
Wait for the user's response. Once they select a phase, provide the exact steps they need to follow based on the `MasterSOP.md`:

**If they choose Phase -1 (Legacy Onboarding):**
1. Run `/discover` to scan the codebase and populate `spec/compiled/architecture.yml`.
2. Review the generated YAML map for accuracy.
3. Run `/baseline` to generate the initial `SuperPRD.md`.
4. Once verified, they can move to Phase 1.

**If they choose Phase 0 (System Initialization):**
1. Ensure the template is cloned with `.agents/` and `.gitkeep` directories.
2. Run `pip install pyyaml`.
3. Give execution permissions to the scripts: `chmod +x .agents/scripts/*.py` (Note: ensure correct path `.agents/scripts/`).

**If they choose Phase 1 (Specification):**
1. Run `/architect` to begin the requirements extraction interview (creates `Draft_PRD.md`).
2. **Open a New Context Window** and run `/redteam` for adversarial analysis (creates `RedTeam_Report.md`).
3. **Open a New Context Window** and run `/resolve` to resolve trade-offs and compile the final `SuperPRD.md` and `MiniPRD` files.
4. Run the memory flush script: `python .agents/scripts/archive_specs.py [Feature_Name]`.

**If they choose Phase 2 (Execution):**
1. **Open a New Context Window** and instruct the Builder agent to implement a specific MiniPRD.
2. Run `python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [modified_node_ids]` after code is written.
3. **Open a New Context Window** and strictly run `/audit spec/compiled/MiniPRD_[Target].md` to verify the code against the contract.

### Step 4: Iterative & Testing Notes
Conclude by reminding the user of two important framework rules:
- **Testing**: Subjective/AI-generated tests go to `tests/candidate_outputs/` for human-in-the-loop review before moving to `tests/fixtures/`.
- **Iterative Updates**: Future features follow the identical Phase 1 ➞ Phase 2 loop, relying on the agents to detect existings YAMLs and perform "Delta Extractions."

### Step 5: Next Actions
Ask the user if they are ready to begin their chosen phase or if they have any specific questions about the slash commands (like `/architect` or `/redteam`).
