# META.yml Template: Skill Model Assignment

**Purpose:** Store model routing metadata for skills to enable intelligent model tier selection.

---

## Template Structure

```yaml
---
name: [skill_name]
description: "[One-line description of skill function]"

# Model assignment: haiku | sonnet | opus
assigned_model: [model_tier]

# Optional: pin to specific model version
# If omitted, uses current version from heuristic_rules.yaml
# model_version: "claude-haiku-4-5-20251001"

# Optional: override max thinking tokens
# max_thinking_tokens: 8192

# Optional: metadata for heuristic classifier (helps auto-suggest tiers)
# tool_calls: 5
# reasoning_intensity: "medium"  # none | low | medium | high | very_high
# output_determinism: "high"  # high | medium | low

# Optional: cost/latency indicators
# cost_sensitivity: "low"  # low | medium | high
# latency_critical: false
```

---

## Field Descriptions

| Field | Required | Type | Description |
|---|---|---|---|
| **name** | Yes | string | Skill identifier (e.g., "hyper-architect") |
| **description** | Yes | string | One-line description of skill purpose |
| **assigned_model** | Yes | enum | Model tier: haiku, sonnet, opus |
| **model_version** | No | string | Specific model ID (e.g., "claude-opus-4-7") |
| **max_thinking_tokens** | No | integer | Override default thinking token limit |
| **tool_calls** | No | integer | Expected number of tool invocations |
| **reasoning_intensity** | No | enum | Reasoning depth: none, low, medium, high, very_high |
| **output_determinism** | No | enum | Output consistency: high, medium, low |
| **cost_sensitivity** | No | enum | Cost vs. quality preference: low, medium, high |
| **latency_critical** | No | boolean | Whether latency is a constraint |

---

## Examples

### Example 1: Read-Only Analysis (Haiku)

```yaml
---
name: hyper-status
description: Display current Living Master Plan and project status

assigned_model: haiku
model_version: "claude-haiku-4-5-20251001"
max_thinking_tokens: 1024

tool_calls: 2
reasoning_intensity: "low"
output_determinism: "high"
```

**Rationale:** Status display is deterministic, minimal reasoning, 2-3 file reads. Haiku is fast and cost-effective.

---

### Example 2: Code Analysis (Sonnet)

```yaml
---
name: hyper-audit
description: Verify written code against its strict requirements and reconcile the Hypergraph memory

assigned_model: sonnet
model_version: "claude-sonnet-4-6"
max_thinking_tokens: 8192

tool_calls: 5
reasoning_intensity: "medium"
output_determinism: "medium"
cost_sensitivity: "medium"
```

**Rationale:** Requires moderate reasoning, multiple file reads and edits, YAML reconciliation. Sonnet balances speed and accuracy.

---

### Example 3: Adversarial Analysis (Opus)

```yaml
---
name: hyper-redteam
description: Stress-test the Draft PRD for security, scalability, and logic flaws

assigned_model: opus
model_version: "claude-opus-4-7"
max_thinking_tokens: 100000

tool_calls: 10
reasoning_intensity: "very_high"
output_determinism: "low"
cost_sensitivity: "low"
latency_critical: false
```

**Rationale:** Adversarial analysis requires deep reasoning and novel problem-solving. Opus is necessary for finding non-obvious vulnerabilities. Cost and latency are not constraints.

---

### Example 4: Generic Skill (Auto-Classified)

```yaml
---
name: hyper-new-skill
description: Example skill awaiting classification

# assigned_model: (will be suggested by heuristic classifier)
```

**Workflow:**
1. Developer creates skill with description only
2. `heuristic_classifier.py suggest()` analyzes description and content
3. Suggests appropriate tier with confidence score
4. Developer confirms or overrides
5. Classifier writes suggested tier to META.yml

---

## Routing Behavior

### No META.yml
If a skill has no META.yml file:
1. ModelRouter detects missing file
2. Logs warning: `{ level: "WARN", event: "opus_fallback", skill: "...", reason: "missing_metadata" }`
3. Emits terminal alert: `⚠️ [skill] has no model assignment. Routing to Opus (fallback).`
4. Routes to Opus (safest option for unknown complexity)

### User Override
User can temporarily override assignment for a single run:
```bash
/hyper-config set-model /hyper-execute sonnet --scope single_run
```

Result: Executes on Sonnet; reverts to META.yml assignment after execution.

---

## Version Pinning

Meta.yml can pin to a specific model version:

```yaml
assigned_model: haiku
model_version: "claude-haiku-4-5-20251001"
```

ModelRouter checks version compatibility:
- ✅ Supported version → Routes normally
- ⚠️ Retired version → Logs warning, routes with alert
- ❌ Invalid version → Falls back to current version

---

## Auto-Classification Workflow

For new skills without model assignment:

```python
classifier = HeuristicClassifier()
spec = SkillSpec(
    name="hyper-new-skill",
    description="My new skill description"
)
suggestion = classifier.suggest(spec)

# Output:
# Suggested Tier: sonnet
# Confidence: 0.88
# Reasoning:
#   - Tool calls: 5
#   - Reasoning intensity: medium
#   - Output determinism: medium
#   - Matched rule: multi_tool_coordination
```

Developer can:
1. **Accept**: Write `assigned_model: sonnet` to META.yml
2. **Override**: Write `assigned_model: opus` instead
3. **Defer**: Leave empty; system falls back to Opus

---

## Cost/Performance Trade-offs

Model tiers by cost and latency:

| Tier | Cost/MTok | Latency | Max Tokens | Use When |
|---|---|---|---|---|
| **Haiku** | $0.80 | ~100ms | 8K | Read-only, simple text |
| **Sonnet** | $3.00 | ~500ms | 8K | Code gen, moderate reasoning |
| **Opus** | $15.00 | ~2000ms | 200K | Complex reasoning, novel problems |

Assignment decision tree:
- Simple, deterministic work? → **Haiku**
- Moderate complexity, code generation? → **Sonnet**
- Complex reasoning, adversarial analysis? → **Opus**
- Unknown? → **Opus** (fallback; safest choice)

---

## Integration Points

1. **Skill Creation**: Create META.yml alongside skill definition
2. **/hyper-execute Start**: ModelRouter reads META.yml, routes to appropriate model
3. **User Override**: `/hyper-config set-model` temporarily changes assignment
4. **Version Check**: ModelRouter warns if pinned version is retired
5. **Auto-Classification**: HeuristicClassifier suggests tier for new skills
