---
name: hyper-clear
description: Flushes conversation history and resets session context between feature cycles. Idempotent—safe to call multiple times.
---

# Clear Context

This skill clears session context to prevent cross-feature bleed and allow fresh feature cycles. It is idempotent and safe to invoke at any time.

## When to use this skill

- **Automatically** after `/hyper-audit` completes (post-audit hook).
- **Manually** at any time to reset session context between feature cycles.
- Whenever you want to guarantee a fresh session without old conversation history contaminating the next phase.

---

## How to use it

### Step 1: Check Idempotency State

Read `.agents/memory/activeContext.md`:
- If the `session_cleared` flag is `true`: /hyper-clear has already been called in this session. Return OK (idempotent, no action needed).
- If the `session_cleared` flag is `false` or missing: Proceed to Step 2.

### Step 2: Flush Conversation History

Clear conversation history by **starting a new conversation** in Claude Code:
1. Inform the user: "Context cleared. Starting fresh session for next feature cycle."
2. The user manually starts a new conversation thread (or you can suggest: "You may now `/hyper-clear` by opening a new chat window").

**Specification Retained:**
- `spec/compiled/architecture.yml` — preserved
- `spec/compiled/SuperPRD.md` — preserved
- `spec/compiled/MiniPRD_*.md` — preserved
- `.agents/` framework code — preserved
- Performance metrics (heuristics, model overrides) — preserved

**Session State Flushed:**
- Conversation history from this context window — cleared
- Temporary working memory — cleared

### Step 3: Update Idempotency Tracker

Edit `.agents/memory/activeContext.md`:
```markdown
## Session State
- session_cleared: true
- cleared_at: [CURRENT_TIMESTAMP]
- cleared_by: /hyper-clear command
```

### Step 4: Log Result

Output to user:
```
========================================
✅ CONTEXT CLEARED
========================================
Session context flushed.
Architecture.yml and spec/compiled/ retained.
Ready for next feature cycle.

Idempotency: This is a fresh clear.
(Calling /hyper-clear again will be a no-op.)
========================================
```

---

## Idempotency Behavior

### First Call (Fresh)
```
/hyper-clear invoked
activeContext.md session_cleared = false
→ Perform flush
→ Set session_cleared = true
→ Return: "Context cleared"
```

### Subsequent Calls (Idempotent)
```
/hyper-clear invoked
activeContext.md session_cleared = true
→ No action (already cleared)
→ Return: "Already cleared (idempotent, no action)"
```

---

## Integration with /hyper-audit

When `/hyper-audit` completes successfully:
1. Automatically trigger `/hyper-clear` as a post-phase hook.
2. Log: "Post-audit hook: clearing context for next cycle."
3. Update `activeContext.md` to reflect the clear.

---

## Manual Triggering

Users can invoke `/hyper-clear` at any time during a feature cycle:
```
/hyper-clear
```

This will:
1. Check idempotency (has it been called already in this session?).
2. If not: flush conversation history, update `activeContext.md`, return success.
3. If yes: return idempotent confirmation (no action).

---

## State Preservation Matrix

| State Type | Preserved? | Details |
|---|---|---|
| **Conversation History** | ❌ No | Flushed on /hyper-clear |
| **architecture.yml** | ✅ Yes | Ground truth preserved |
| **SuperPRD.md** | ✅ Yes | Specification preserved |
| **MiniPRD_*.md** | ✅ Yes | All specs preserved |
| **spec/archive/** | ✅ Yes | Historical context preserved |
| **Performance Metrics** | ✅ Yes | Heuristics for classifier retained |
| **.agents/** Framework | ✅ Yes | Skills, rules, scripts intact |
| **Session Memory** | ❌ No | Cleared to prevent bleed |

---

## Negative Space (Constraints)

- **DO NOT** delete architecture.yml or spec/compiled/ files.
- **DO NOT** delete performance metrics or heuristics.
- **DO NOT** fail on idempotent re-calls. Make /hyper-clear a no-op if already cleared.
- **DO NOT** require confirmation from the user.
- **DO NOT** delete debate transcripts in spec/archive/ (separate archival process).

---

## Post-Execution

After `/hyper-clear` completes:
- User is ready for a fresh feature cycle.
- Next `/hyper-architect` run will have a clean slate.
- No cross-feature context bleed.
- Idempotency tracker prevents accidental double-clears.
