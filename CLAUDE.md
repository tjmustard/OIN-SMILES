# Claude Code Integration Guide

> **For Human Developers:** See `README.md` for the complete usage guide, available slash
> commands, and setup instructions.

---

## For Claude Code: System Mandates

Read `AGENTS.md` for the full framework system mandates that apply to all IDEs.

**Claude Code-specific overrides and additions:**

### Tool Names
When skills say "read/write/run/edit a file," use these Claude Code tools:

| Action | Tool |
|---|---|
| Read a file | **Read** tool |
| Write a file | **Write** tool |
| Edit a file | **Edit** tool |
| Run a shell command | **Bash** tool |
| Ask the user a question | **AskUserQuestion** tool |
| Search file patterns | **Glob** tool |
| Search file contents | **Grep** tool |

### Skill Invocation
Skills are available as slash commands in `.claude/commands/`. Each command is a thin bridge
that reads its corresponding `.agents/skills/hyper-<name>/SKILL.md`. When the user invokes a command,
the skill file provides the full instructions.

### Task Tracking
For any task involving 3 or more steps, use the built-in task tools **before** starting work:
- **TaskCreate** — create tasks with clear subjects
- **TaskUpdate** — mark `in_progress` when starting, `completed` when done
- **TaskList** — check overall progress

### Context Window Management
When a skill instructs you to "open a new context window": **complete the current agent turn**,
then inform the user to start a new conversation thread for the next phase. This prevents
cross-contamination between adversarial agents (the Red Team must not see the Architect's
conversation history).

---

## Schema Definitions (Persistent Rules)

**Note:** These schemas are embedded in CLAUDE.md to reduce per-message token overhead.
For the authoritative versions, see `.agents/schemas/` directory (deprecated; migration in progress).

### SuperPRD Schema

SuperPRD (Super Product Requirements Document) is the comprehensive specification for a feature.
Use this template when compiling a full feature specification from architect and red team input.

**Structure:**
```
# SuperPRD: [Feature Name]

## 1. Introduction & Goals
- Problem Statement: [Why are we building this?]
- Solution Overview: [High-level approach]
- Target Audience: [Who is this for?]

## 2. Confidence Mandate
- Confidence Score: [1-10] (must be calculated before proceeding)
- Clarifying Questions: [List any open questions]

## 3. Scope
- In-Scope: [Features included]
- Out-of-Scope: [Explicitly excluded features]

## 4. User Stories (Atomic)
| ID | User Story | Acceptance Criteria | Priority |
| US-001 | As [User], I want [Action] so that [Outcome] | 1. Criterion A<br>2. Criterion B | High |

## 5. Technical Specifications
- Architecture & Resolved Trade-offs: [System design + trade-off log]
- System Graph Blast Radius: [Affected nodes in architecture.yml]
- Execution Checklist: [List of MiniPRDs to execute]
- API Contracts / Schema: [Type definitions]
- Dependencies: [Libraries/frameworks]

## 6. Negative Constraints
- **DO NOT** [Anti-pattern 1]
- **DO NOT** [Anti-pattern 2]

## 7. Risks & Mitigation
- Risk 1: [Description] → Mitigation: [Action]

## 8. Success Metrics
- [Metric 1]
- [Metric 2]
```

**Reference:** See `.agents/schemas/SuperPRD_Template.md` (deprecated; use CLAUDE.md above)

---

### MiniPRD Schema

MiniPRD (Mini Product Requirements Document) is a modular, executable specification for a single
feature module. Generate one MiniPRD per independent task.

**Structure:**
```
# MiniPRD: [Module Name]

**Hypergraph Node ID:** [node_id]
**Parent Node:** [parent_node_id]

## 1. The Confidence Mandate
- Confidence Score: [1-10] (required before implementation)
- Clarifying Questions: [If < 9, list questions needed]

## 2. Atomic User Stories
- US-001: As [User Type], I want [Action] so that [Value]
- US-002: ...

## 3. Implementation Plan (Task List)
- [ ] Task 1: [Specific, <10 min effort]
- [ ] Task 2: ...

## 4. The Negative Space (Constraints)
- **DO NOT** [Anti-pattern]
- **DO NOT** [Architectural violation]

## 5. Integration Tests & Verification
- Test 1 (Deterministic): [Input] → [Expected Output]
- Test 2 (Novel): [Input] → [Candidate Artifact routing]
```

**Reference:** See `.agents/schemas/MiniPRD_template.md` (deprecated; use CLAUDE.md above)

---

### Hypergraph Schema (architecture.yml)

The hypergraph is a YAML file (`spec/compiled/architecture.yml`) that tracks system dependencies as a directed acyclic graph.

**Node Structure:**
```yaml
nodes:
  - id: [unique_identifier]
    dimension: [System | Module | Atomic]  # Layer of abstraction
    status: [clean | dirty | needs_review]  # Build state
    associated_file: [path_to_source]  # MiniPRD, source code, or doc
    description: [semantic_purpose]  # What this node does
    inputs:
      - data_type: [type_name]
        source_id: [upstream_node_id]
    outputs:
      - data_type: [type_name]
        target_id: [downstream_node_id]
    edges:
      depends_on: [list_of_node_ids]  # Architectural dependencies
      implements: [list_of_node_ids]  # Hierarchical link (Atomic→Module)
```

**Status Values:**
- `clean` — Implementation verified against specification; ready for use
- `dirty` — Recently modified; awaiting audit review
- `needs_review` — Dependent on modified node; blast-radius mark

**Reference:** See `.agents/schemas/hypergraph_schema.md` (deprecated; use CLAUDE.md above)

---

## Migration Notice

**Effective:** 2026-05-03  
**Deadline:** 2026-06-17 (6-week deprecation period)

All `.agents/schemas/` files have been migrated to CLAUDE.md to reduce per-message context overhead.
Old files marked with `.deprecated` suffix; will be deleted after 2026-06-17.

**For Developers:**
- Update custom skills to reference "CLAUDE.md: Schema Definitions" instead of `.agents/schemas/` files
- See `docs/MIGRATION_RULES.md` for detailed migration guide
- CI lint will enforce CLAUDE.md-only references post-deadline
