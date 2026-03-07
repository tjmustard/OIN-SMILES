
# SuperPRD Template

## Metadata
- **Project Name**: [Project Name]
- **Version**: 1.0.0
- **Status**: [Draft / Review / Approved]
- **Owner**: [Your Name / Agent Name]

## 1. Introduction & Goals
### 1.1 Problem Statement
[Describe the problem clearly. Why are we building this?]

### 1.2 Solution Overview
[High-level solution description.]

### 1.3 Target Audience
[Who is this for?]

## 2. Confidence Mandate
**Confidence Score**: [1-10] (Initial score must be calculated before proceeding)
**Clarifying Questions**:
- [ ] Question 1
- [ ] Question 2

## 3. Scope
### 3.1 In-Scope
- Feature 1
- Feature 2

### 3.2 Out-of-Scope
- Feature X
- Feature Y

## 4. User Stories (Atomic)
| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As a [User], I want to [Action], so that [Outcome]. | 1. Criterion A<br>2. Criterion B | High |
| US-002 | ... | ... | Medium |

## 5. Technical Specifications (The Blueprint)
### 5.1 Architecture & Resolved Trade-offs
[Describe the system architecture, data flow, etc.]

**Resolved Trade-offs Log:**
- **Issue:** [Description of architectural conflict raised by Red Team]
- **Options Considered:** [Option A] vs [Option B]
- **Resolution:** [Why the specific decision was made]

### 5.2 System Graph Blast Radius
The following nodes in `spec/compiled/architecture.yml` will be affected/modified by this feature:
- `[node_id_1]`
- `[node_id_2]`

### 5.3 Execution Checklist (MiniPRDs)
The following modular PRDs have been generated to execute this SuperPRD in isolated branches:
- [ ] `spec/compiled/MiniPRD_[Module_Name_1].md`
- [ ] `spec/compiled/MiniPRD_[Module_Name_2].md`

### 5.4 API Contracts / Schema
```typescript
interface User {
  id: string;
  email: string;
  // ...
}
```

### 5.5 Dependencies
- [List major libraries/frameworks]

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** use [Library X].
- **DO NOT** expose [Sensitive Data].
- **DO NOT** modify [Existing System Y].

## 7. Risks & Mitigation
- **Risk 1**: [Description] -> **Mitigation**: [Action]

## 8. Success Metrics
- [Metric 1]
- [Metric 2]
