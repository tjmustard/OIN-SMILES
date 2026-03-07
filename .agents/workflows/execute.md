---
description: Execute the implementation plan
---

# /execute Workflow

Now implement precisely as planned, in full.

Implementation Requirements:

- Write elegant, minimal, modular code.
- Adhere strictly to existing code patterns, conventions, and best practices.
- Include thorough, clear comments/documentation within the code.
- As you implement each step:
  - Update the markdown tracking document with emoji status and overall progress percentage dynamically.
- **MANDATORY**: Once you have modified the code, execute `python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [modified_node_ids]` to deterministically flag the hypergraph. Failure to do this corrupts the state-machine.