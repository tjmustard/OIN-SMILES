# Security & Safety

## 1. Input Sanitization (Critical)
- **Data Validation**:
  - ALL incoming data must be validated (e.g. using Pydantic).
- **Untrusted Input**:
  - Validate file formats strictly.
  - Do not execute arbitrary code embedded in input files (pickle, etc.).

## 2. Execution Constraints
- **File System**:
  - Read/Write ONLY within the workspace or explicitly provided output directories.
- **Secrets**:
  - No API keys in source code. Use `.env` or secret management.

## 3. Library Safety
- **External Calls**:
  - Handle external library failures gracefully.
  - Wrap complex operations in try/except blocks to prevent crashes.
