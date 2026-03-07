# Claude Code Integration Guide

This document serves as both a manual for human developers and **foundational system instructions** for Claude Code. By reading this file, Claude Code understands how to natively orchestrate the Hypergraph Coding Agent Framework.

## For Human Developers

The Hypergraph Framework was originally designed for manual LLM copy-pasting and manual script execution. With **Claude Code**, the entire process is automated using built-in tools.

### The Automated Workflow
You do not need to manually open new context windows or run the Python scripts yourself. Simply instruct Claude Code to transition through the phases:
1. **Onboarding**: "Run the Discover phase" or "Run the Baseline phase."
2. **Specification**: "Activate the Architect skill and let's design a new feature." → "Now run the Red Team skill on that draft." → "Let's Resolve the trade-offs."
3. **Execution**: "Implement MiniPRD-1." Claude Code will automatically check the blast radius, write the code, and verify it.

### Available Slash Commands
Use these commands to invoke the framework's specialized agents:

| Command | Phase | Description |
|---|---|---|
| `/architect` | Phase 1 | State-machine interview → generates `Draft_PRD.md` |
| `/redteam` | Phase 1 | Adversarial analysis → generates `RedTeam_Report.md` |
| `/resolve` | Phase 1 | Mediates trade-offs → compiles `SuperPRD.md` + `MiniPRD` files |
| `/audit` | Phase 2 | Verifies code against MiniPRD → reconciles `architecture.yml` |
| `/discover` | Phase -1 | Scans codebase → initializes `architecture.yml` |
| `/baseline` | Phase -1 | Reverse-engineers system → generates baseline `SuperPRD.md` |
| `/sop` | Any | Explains the Master SOP and guides you to the right phase |

---

## For Claude Code: System Mandates

**CRITICAL INSTRUCTIONS FOR CLAUDE CODE:** When operating in this workspace, you must adhere to the following framework rules:

### 1. Skill Integration
You have native access to the framework's phases via the custom slash commands defined in `.claude/commands/`. When the user requests a phase or invokes a command, load and follow its specific procedural guidance strictly.

### 2. Autonomous Tool Execution
- You MUST autonomously execute the deterministic Python scripts in `.agents/scripts/` using your **Bash tool** when mandated by a skill's instructions.
- If a skill requires calculating the blast radius, run:
  ```bash
  python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [node_ids]
  ```
- If a skill requires flushing the active workspace, run:
  ```bash
  python .agents/scripts/archive_specs.py [Feature_Name]
  ```
- Always verify the script completed successfully by checking the exit code and output.

### 3. State Management
- Treat `spec/compiled/architecture.yml` as the **absolute ground truth** for the project's state.
- Write temporary drafts and reports to `spec/active/` using the **Write** tool.
- Ensure all generated specifications strictly adhere to the templates in `.agents/schemas/`.
- Use the **Read** tool to inspect existing specs before modifying or generating new ones.

### 4. Interactive Interviews
When acting as the Architect or Resolution Agent, use the **AskUserQuestion** tool to gather structured constraints efficiently. Enforce the **Pacing Loop**: never ask more than **2 questions per turn**.

### 5. No Probabilistic Traversal
Do not guess architectural dependencies. Always rely on the output of `hypergraph_updater.py` to understand the blast radius before executing code modifications.

### 6. File Lifecycle Rules
- **Read** existing files before overwriting them.
- **Never** write to `spec/archive/` manually — use `archive_specs.py` exclusively.
- **Never** read from `spec/archive/` or `tests/candidate_outputs/` during agentic tasks (treat these as blocked per `.agentignore`).

### 7. Rules Always Active
Apply the project rules from `.agents/rules/` to all code generation:
- `.agents/rules/python.md` — Python style, type hints, and architectural constraints
- `.agents/rules/security.md` — Input validation, secrets management, file system constraints
- `.agents/rules/testing.md` — Testing standards
- `.agents/rules/package-management.md` — Dependency management

### 8. Task Tracking
For any task involving 3 or more steps, always use the built-in task list tools (TaskCreate, TaskUpdate, TaskList) to track progress before starting work.

### 9. Context Window Management
When a skill instructs you to "open a new context window," in Claude Code this means: **complete the current agent turn**, then inform the user to start a new conversation thread for the next phase. This prevents cross-contamination between adversarial agents (e.g., the Red Team must not see the Architect's conversation).

---

## Directory Structure Reference

```
.agents/
├── skills/         # Skill definitions (source of truth for .claude/commands/)
├── schemas/        # Immutable templates (MiniPRD, SuperPRD, hypergraph)
├── scripts/        # Deterministic state management
│   ├── hypergraph_updater.py   # BFS blast radius propagation
│   └── archive_specs.py        # Active spec flushing
├── workflows/      # Workflow definitions
├── rules/          # Always-on coding standards
└── memory/         # Project context files

spec/
├── active/         # Working drafts (temporary, will be archived)
├── compiled/       # Ground truth (SuperPRD, MiniPRDs, architecture.yml)
└── archive/        # Historical — DO NOT READ during agent tasks

tests/
├── candidate_outputs/  # Unverified AI outputs — DO NOT READ during tasks
└── fixtures/           # Verified regression baselines
```
