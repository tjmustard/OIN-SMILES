# Final Truth Template: Context Compaction Output

**Purpose:** Preserve decision reasoning from `/hyper-resolve` while discarding verbatim debate.

---

## Template Structure

```markdown
# Final Truth: [Feature Name]

**Generated:** [ISO 8601 timestamp]
**Feature:** [Feature name from /hyper-resolve]
**Archive:** See full transcript: `spec/archive/[YYYY-MM-DD]_[8-char-hash]_debate.md`

---

## Trade-off Decisions

### 1. [Trade-off Title]

**Architect Position:**
[Concise statement of architect's preferred approach + rationale]

**Red Team Position:**
[Concise statement of red team's counter-argument + why]

**Resolution:**
[Final decision made + reasoning for why this path was chosen]

**Risks Identified:**
- [Risk 1]
- [Risk 2]
- ...

**Mitigations:**
- [Mitigation 1]
- [Mitigation 2]
- ...

---

### 2. [Next Trade-off Title]
[Repeat structure above...]

---

## Archive & Retention

Full debate transcript archived at: `spec/archive/[YYYY-MM-DD]_[hash]_debate.md`

**Retention Policy:** Minimum 90 days. This transcript is the authoritative record
of all decisions made and rejected during this feature cycle.

**Link back to compiled specs:** See `spec/compiled/` for SuperPRD and MiniPRDs.
```

---

## Field Descriptions

| Field | Required | Description |
|---|---|---|
| **Feature Name** | Yes | Name of the feature being resolved |
| **Generated Timestamp** | Yes | ISO 8601 UTC timestamp of generation |
| **Archive Path** | Yes | Reference to full debate transcript (for auditing) |
| **Trade-off Title** | Yes | Concise 2-4 word title of the decision |
| **Architect Position** | Yes | Architect's preferred approach (1-3 sentences) |
| **Red Team Position** | Yes | Red Team's counter-argument (1-3 sentences) |
| **Resolution** | Yes | Final decision and reasoning (2-4 sentences) |
| **Risks** | Conditional | Include if trade-off involved risk assessment |
| **Mitigations** | Conditional | Include if mitigations were identified |
| **Retention Note** | Yes | Always include 90-day minimum retention policy |

---

## Constraints (Negative Space)

- **DO NOT** include verbatim debate text (that's in the archive)
- **DO NOT** include cancelled ideas or rejected options (only what was decided)
- **DO NOT** include preliminary explorations or brainstorming
- **DO NOT** omit the Architect and Red Team positions (both are required context)
- **DO NOT** omit the Archive path reference (needed for auditing decisions)

---

## Example: Filled Template

```markdown
# Final Truth: Context Compaction

**Generated:** 2026-05-03T12:30:45Z
**Feature:** Context Compaction & Trade-off Preservation
**Archive:** See full transcript: `spec/archive/2026-05-03_abc123de_debate.md`

---

## Trade-off Decisions

### 1. Compaction Aggressiveness

**Architect Position:**
Aggressive compaction: extract only final decisions, discard all debate. Minimizes context
for next team.

**Red Team Position:**
Conservative compaction: preserve both positions + resolution. Developers need to
understand "why" decisions were made, not just "what" was decided.

**Resolution:**
Conservative approach chosen. Preserve Architect and Red Team positions with their
reasoning. Archive full transcript for auditing. This prevents "mysterious decisions"
downstream.

**Risks Identified:**
- Final Truth document could still be large if many trade-offs exist
- Developers might miss important nuances by reading condensed version

**Mitigations:**
- Final Truth includes link to full transcript for deeper reading
- MiniPRD design should reference Final Truth for context

---

### 2. Archive Retention Duration

**Architect Position:**
Keep archives indefinitely. They're small, valuable for historical auditing.

**Red Team Position:**
90-day retention minimum. Older archives crowd spec/archive/, create confusion,
require maintenance.

**Resolution:**
90-day minimum retention with cleanup script. Older archives can be moved to
cold storage if needed, but spec/archive/ stays manageable.

**Risks Identified:**
- 90 days might be too short for compliance auditing in some domains
- Cleanup script needs to be reliable

**Mitigations:**
- Policy can be updated per-project in .agents/scripts/interrupt_detector.py
- Cleanup only runs on user request or scheduled maintenance

---

## Archive & Retention

Full debate transcript archived at: `spec/archive/2026-05-03_abc123de_debate.md`

**Retention Policy:** Minimum 90 days. This transcript is the authoritative record
of all decisions made and rejected during this feature cycle.

**Link back to compiled specs:** See `spec/compiled/` for SuperPRD and MiniPRDs.
```

---

## Output File Naming

Format: `spec/compiled/final_truth_[YYYY-MM-DD]_[Feature_Name].md`

Example: `spec/compiled/final_truth_2026-05-03_ContextCompaction.md`

---

## Integration Points

1. **Post-/hyper-resolve Hook:** Automatically generate Final Truth after `/hyper-resolve` completes
2. **Pre-/hyper-execute Check:** `interrupt_detector.py` verifies Final Truth exists; auto-recover if missing
3. **Archival:** Full transcript moved to `spec/archive/` with metadata and retention reminder
4. **MiniPRD Generation:** MiniPRDs should reference Final Truth for context on key trade-offs
