# Hypergraph Coding Agent Framework — Agent Instructions

> **HACF as a Toolchain:** This project uses the Hypergraph Coding Agent Framework
> (HACF) as its development toolchain. The skills in `.agents/skills/`, the scripts
> in `.agents/scripts/`, and the schemas in `.agents/schemas/` are development tools —
> they are **not** subjects of this project's plans, PRDs, or architecture docs.
> When you create SuperPRDs, MiniPRDs, or architecture nodes, you are documenting
> **this project**, not the HACF framework itself.
>
> Cross-IDE manifest: injected as always-on context by Windsurf, Cursor, Roo Code,
> GitHub Copilot, and Zed. IDE-specific bridge files in `.claude/`, `.cursor/`,
> `.windsurf/`, `.clinerules/`, and `.roo/` extend this base with tool-specific overrides.
> For the human-facing usage guide, see `README.md`.

---

## Single Source of Truth

`.agents/` is the authoritative location for **all** skill, rule, schema, and script content.
All other IDE directories contain only thin bridge files that reference `.agents/`.

```
.agents/
├── skills/     # All skill/command definitions — one directory per skill
├── schemas/    # Immutable templates (MiniPRD, SuperPRD, hypergraph)
├── scripts/    # Deterministic state tools (hypergraph_updater.py, archive_specs.py)
├── rules/      # Always-on coding standards
└── memory/     # Project context files (activeContext, productContext, systemPatterns)

spec/
├── active/     # Working drafts — TEMPORARY, will be archived
├── compiled/   # Ground truth (SuperPRD, MiniPRDs, architecture.yml)
└── archive/    # Historical — BLOCKED from agent reads (.agentignore)

tests/
├── candidate_outputs/  # Unverified AI outputs — BLOCKED from agent reads
└── fixtures/           # Verified regression baselines
```

---

## System Mandates

### 1. Skill Invocation
When the user invokes a skill command (e.g., `/hyper-architect`, `/hyper-redteam`), read the
corresponding `.agents/skills/<name>/SKILL.md` and follow its instructions precisely.

To discover all available skills:
```
ls .agents/skills/
```

### 2. Autonomous Script Execution
Execute these deterministic Python scripts when mandated by a skill:

```bash
# Propagate blast radius after modifying code:
python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [node_id_1] [node_id_2]

# Flush active workspace after specification phase completes:
python .agents/scripts/archive_specs.py [Feature_Name]
```

Always verify the script completed successfully by checking the exit code and output.

### 3. State Management
- `spec/compiled/architecture.yml` is the **absolute ground truth** for project state.
- Write working drafts and reports to `spec/active/`.
- All generated specs must strictly follow the templates in `.agents/schemas/`.
- Read existing files before overwriting them.

### 4. Interactive Interviews
When acting as the Architect or Resolution Agent, enforce the **Pacing Loop**:
ask a MAXIMUM of **2 questions per turn**. Wait for the user's response before proceeding.

### 5. No Probabilistic Traversal
Never guess architectural dependencies. Always rely on `hypergraph_updater.py` output
to understand the blast radius before executing code modifications.

### 6. File Lifecycle Rules
- **Never** write to `spec/archive/` manually — use `archive_specs.py` exclusively.
- **Never** read from `spec/archive/` or `tests/candidate_outputs/` during agentic tasks
  (treated as blocked per `.agentignore`).
- **Never** edit `spec/compiled/architecture.yml` directly — use `/hyper-audit` or `/hyper-discover`.

### 7. Always-On Coding Rules
Apply the rules in `.agents/rules/` to all code generation:

| Rule File | Scope |
|---|---|
| `python.md` | Python style, type hints, architectural constraints |
| `security.md` | Input validation, secrets management, file system constraints |
| `testing.md` | Testing standards and patterns |
| `package-management.md` | Dependency management with `uv` |

---

## Available Skills

| Skill | Trigger | Phase | Description |
|---|---|---|---|
| `hyper-architect` | `/hyper-architect` | 1 | Requirements interview → Draft_PRD.md |
| `hyper-redteam` | `/hyper-redteam` | 1 | Adversarial analysis → RedTeam_Report.md |
| `hyper-resolve` | `/hyper-resolve` | 1 | Trade-off mediation → SuperPRD + MiniPRDs |
| `hyper-audit` | `/hyper-audit` | 2 | Code verification → reconciles architecture.yml |
| `hyper-execute` | `/hyper-execute` | 2 | Checks activeContext.md, implements a MiniPRD, updates hypergraph |
| `hyper-discover` | `/hyper-discover` | -1 | Scans codebase → initializes architecture.yml |
| `hyper-baseline` | `/hyper-baseline` | -1 | Reverse-engineers system → baseline SuperPRD |
| `hyper-init` | `/hyper-init` | 0 | Scaffolds standard repository documentation templates |
| `hyper-clear` | `/hyper-clear` | post-audit | Idempotent context flushing; flushes conversation history while preserving specs and metrics |
| `hyper-contextualize` | `/hyper-contextualize` | any | Audit and fix installed project agent files (`CLAUDE.md`, `AGENTS.md`, `GEMINI.md`) to ensure HACF is framed as a development toolchain, preventing framework content from bleeding into project plans and architecture docs |
| `hyper-sop` | `/hyper-sop` | any | Master SOP guide and phase orientation |
| `hyper-status` | `/hyper-status` | any | Living Master Plan snapshot |
| `hyper-consult-cto` | `/hyper-consult-cto` | pre-spec | CTO advisor for architectural decisions |
| `hyper-co-research` | `/hyper-co-research` | any | Peer-level AI research partner |
| `hyper-deepdive` | `/hyper-deepdive` | any | Exhaustive First Principles topic research |
| `hyper-create-skill` | `/hyper-create-skill` | any | Convert a prompt into a new skill |
| `hyper-new-workflow` | `/hyper-new-workflow` | any | Scaffold a new skill and IDE bridges |
| `hyper-document` | `/hyper-document` | any | Update README, CHANGELOG, docs/, AGENTS.md, and memory files after any change |
| `hyper-grill-docs` | `/hyper-grill-docs` | any | Relentless domain-sharpening interview: challenges terminology against CONTEXT.md, cross-references code, stress-tests domain boundaries with concrete scenarios, and writes ADRs inline when decisions are hard-to-reverse, surprising, and trade-off-driven |
| `hyper-handoff` | `/hyper-handoff` | any | Compact the current conversation into a handoff document (saved to OS temp dir) so a fresh agent can continue the work. Pass an optional argument describing what the next session will focus on. |
| `hyper-session-update` | `/hyper-session-update` | any | Sync memory with session work |
| `hyper-update` | `/hyper-update` | any | Smart upgrade framework while preserving customizations |
| `hyper-refresh-memory` | `/hyper-refresh-memory` | any | Rebuild mental model from codebase |
| `hyper-troubleshooting` | `/hyper-troubleshooting` | any | Diagnose framework failure states |
| `hyper-tutorial` | `/hyper-tutorial` | any | Framework walkthrough for new users |
| `hyper-tutorial-generator` | `/hyper-tutorial-generator` | any | Generate markdown tutorial from integration test or user-provided files via iterative LLM collaboration |
| `hyper-stitch-design` | `/hyper-stitch-design` | any | UI/UX → Design System specification |
| `hyper-prompt-engineer` | `/hyper-prompt-engineer` | any | Collaborative prompt design |
| `hyper-template-architect` | `/hyper-template-architect` | any | Reverse-engineer document into template |
| `hyper-peer-review` | `/hyper-peer-review` | any | Evaluate and triage peer review findings |
| `hyper-create-issue` | `/hyper-create-issue` | any | Format and file a GitHub issue |
| `hyper-publish` | `/hyper-publish` | any | AI-proposes commit messages from CHANGELOG; HITL approval gates at each step; optional PR creation with auto-drafted descriptions. Touches git-changed files (mtime), commits, pushes, and optionally creates GitHub PR. |
| `hyper-learning-opportunity` | `/hyper-learning-opportunity` | any | Structured teaching on any concept |

Full skill instructions: `.agents/skills/hyper-<name>/SKILL.md`

---

## Cost Optimization: Hybrid Model Orchestration

### Model Routing

Skills automatically route to the optimal Claude model tier based on metadata:

- **Haiku** (2k thinking ceiling): Routine, deterministic tasks. ~70% cost savings.
- **Sonnet** (10k thinking ceiling): Tactical reasoning, trade-off analysis. ~50% cost savings.
- **Opus** (20k thinking ceiling): Strategic decisions, adversarial analysis.

**Configuration**: Each skill has a `.agents/skills/hyper-<name>/META.yml` file specifying its assigned model:

```yaml
assigned_model: haiku
model_version: "claude-haiku-4-5-20251001"
max_thinking_tokens: 2000  # Optional per-skill override
```

### Manual Overrides

Temporarily override a skill's model assignment:

```bash
# Single execution only (resets after)
/hyper-config set-model /hyper-execute sonnet --scope single_run

# Permanent (until manually reverted)
/hyper-config set-model /hyper-execute sonnet --scope permanent --reason "Debugging complex case"

# Revert to default
/hyper-config revert-override /hyper-execute
```

### Token Enforcement

- **Thinking Token Ceilings**: Enforced per model; execution halts with warning if exceeded. Per-skill overrides available via META.yml.
- **Output Token Budgets**: 50k ceiling per task. MiniPRDs exceeding this are flagged for splitting.
- **Variance Tracking**: Post-execution reconciliation logs estimated vs. actual tokens. Deviations >20% flagged for investigation.

---

## Framework Workflow (Quick Reference)

```
Phase -1 (Brownfield): /hyper-discover → /hyper-baseline
Phase  1 (Spec):       /hyper-architect → /hyper-redteam → /hyper-resolve → archive_specs.py
Phase  2 (Build):      /hyper-execute → hypergraph_updater.py → /hyper-audit
Phase  3 (Novel):      Human review → tests/fixtures/ → update MiniPRD
```

Each phase boundary requires a **fresh context window** to prevent cross-contamination
between adversarial agents (Red Team must not see Architect's conversation history).

---

## Schema Definitions (Persistent Rules)

**Note:** These schemas define the document structures for **this project's** specifications.
Use them when running `/hyper-architect`, `/hyper-resolve`, and `/hyper-execute`.

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
