# ⚠️ DEPRECATED: Hypergraph Schema

**Status:** Migrated to CLAUDE.md (See CLAUDE.md: Schema Definitions › Hypergraph Schema)  
**Deadline for Removal:** 2026-06-17 (6-week deprecation period)

This file is maintained for backward compatibility only. New references should use CLAUDE.md.

---

# architecture.yml Standard
nodes:
  - id: string # Unique identifier (e.g., "auth_module", "login_function")
    dimension: enum [System | Module | Atomic]
    status: enum [clean | dirty | needs_review] # For deterministic traversal script
    associated_file: string # Path to MiniPRD (Module), source code (Atomic), or systemPatterns.md (System)
    description: string # Semantic purpose
    inputs: 
      - data_type: string
        source_id: string # ID of the node providing this input
    outputs:
      - data_type: string
        target_id: string # ID of the node receiving this output
    edges:
      depends_on: [list of node_ids] # Architectural dependencies (libraries, DBs)
      implements: [list of node_ids] # Hierarchical link (Atomic nodes implement Module nodes)