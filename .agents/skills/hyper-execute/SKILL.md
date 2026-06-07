---
name: execute
description: Implements the current plan precisely and in full, following existing code patterns and updating the hypergraph after every modification. Use when ready to write code against a compiled MiniPRD.
---

# Execute

This skill implements a plan precisely, minimally, and in full — adhering to existing conventions, updating progress tracking, and ensuring the hypergraph is kept in sync after every modification.

## When to use this skill

- When the user has a compiled `MiniPRD_*.md` in `spec/compiled/` and is ready to build.
- When the user explicitly runs `/hyper-execute` or asks to "implement the plan."
- After `/hyper-resolve` has compiled the final MiniPRDs and the workspace has been archived.

## How to use it

1. **Memory Check**
   - Read `.agents/memory/activeContext.md` before anything else.
   - Check whether the target MiniPRD is listed as complete or already audited.
   - **If the MiniPRD is marked complete or audited:** HALT. Display the warning, then use **AskUserQuestion**:

     ```
     ⚠️ [MiniPRD name] appears already complete per activeContext.md.
     Run /hyper-audit to archive it, or move it manually:
       mv spec/compiled/[MiniPRD].md spec/archive/[MiniPRD]_AUDITED.md

     Confirm you want to re-implement it?

     - Option A: Yes, re-implement — proceed even though it appears done
     - Option B: No, stop — I will archive or review it instead
     ```

2. **Read the MiniPRD**
   - Read the target `MiniPRD_*.md` from `spec/compiled/`.
   - Confirm a Confidence Score (1–10) against the MiniPRD's "Confidence Mandate" section. If below 9, list the clarifying questions, then use **AskUserQuestion**:

     ```
     Confidence is [N]/10. How do you want to proceed?

     - Option A: Answer my questions — provide clarifications so I can proceed
     - Option B: Proceed anyway — accept the risk and implement with current information
     - Option C: Stop — halt until the MiniPRD is revised
     ```

2.5. **Codebase Orientation (Haiku Sub-Agent)**
   - Before implementing, launch a read-only sub-agent to pre-load all relevant file contents:
     - Use the Agent tool with `subagent_type: "general-purpose"` and `model: "haiku"`
     - Prompt the sub-agent: "Read `spec/compiled/[MiniPRD filename]` in full. From the Implementation Plan task list, identify every file path that must be created or modified. For each file that already exists in the codebase, read its full content and return it verbatim. For new files, list the path and expected purpose only. Return a structured report mapping file paths to their current content."
   - Wait for the sub-agent result before proceeding to Step 3.
   - Use the returned file contents as your working context for implementation.
   - **Fallback:** If the sub-agent returns an error or empty result, proceed to Step 3 and read files directly as normal.
   - Do not re-read files the sub-agent already returned unless you need to verify a specific line.

3. **Implement Precisely**
   - Write elegant, minimal, modular code.
   - Adhere strictly to existing code patterns, conventions, and best practices observed in the codebase.
   - Include clear comments where logic is non-obvious.
   - As you implement each task in the MiniPRD's task list, update the tracking document with emoji status and an overall progress percentage.

4. **Update the Hypergraph (MANDATORY)**
   - Once code modifications are complete, you MUST execute the hypergraph updater via the shell:
     ```bash
     python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [modified_node_ids]
     ```
   - Replace `[modified_node_ids]` with the actual node IDs from `architecture.yml` corresponding to the files you changed.
   - **Failure to do this corrupts the state-machine.** The subsequent `/hyper-audit` agent depends on these flags.

5. **Halt and Report**
   - Output a summary of all files modified and nodes flagged.
   - Instruct the user: "Implementation complete. Start a new conversation and run `/hyper-audit spec/compiled/MiniPRD_[Target].md` to verify the code against the contract."
   - Note: a successful `/hyper-audit` run will move the MiniPRD to `spec/archive/` automatically, removing it from future execute runs.