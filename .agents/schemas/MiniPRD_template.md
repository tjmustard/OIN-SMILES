# ⚠️ DEPRECATED: MiniPRD Template

**Status:** Migrated to CLAUDE.md (See CLAUDE.md: Schema Definitions › MiniPRD Schema)  
**Deadline for Removal:** 2026-06-17 (6-week deprecation period)

This file is maintained for backward compatibility only. New references should use CLAUDE.md.

---

# MiniPRD: [Module Name]
**Hypergraph Node ID:** [node_id]
**Parent Node:** [System level node_id]

## 1. The Confidence Mandate
**Agent Instruction:** Before generating any plans or writing code, analyze this document and output a Confidence Score (1-10). If the score is below 9, list strictly the clarifying questions needed to reach 10.

## 2. Atomic User Stories
* **US-001:** As a [User Type], I want to [Action] so that [Value/Result].
* **US-002:** ...

## 3. Implementation Plan (Task List)
- [ ] Task 1: [Specific, actionable step taking <10 minutes]
- [ ] Task 2: ...

## 4. The Negative Space (Constraints)
* **DO NOT** [Anti-pattern or deprecated method].
* **DO NOT** [Architectural violation].

## 5. Integration Tests & Verification
* **Test 1 (Deterministic):** [Input] -> Expected Output: [Exact Output]
* **Test 2 (Novel):** [Input] -> Expected Output: [Candidate Artifact routing protocol triggered]