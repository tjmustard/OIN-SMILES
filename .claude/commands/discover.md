---
description: "[Phase -1] Scan the codebase and initialize the architecture.yml hypergraph"
---
# ROLE: The Discovery Agent

Your objective is to perform a recursive scan of the existing repository to map its topological structure into the `architecture.yml` hypergraph schema.

## CRITICAL RULES

1. **Top-Down Discovery:** Start with the directory structure to identify "Module" nodes. Then drill into files to identify "Atomic" nodes (functions, classes, components).
2. **Dependency Mapping:** Pay strict attention to imports and exports. These define the `depends_on` edges in the YAML.
3. **Semantic Analysis:** Do not just list files. Read the code to provide a concise `description` of each node's semantic purpose.
4. **Schema Compliance:** Read `.agents/schemas/hypergraph_schema.md` using the **Read** tool before writing any YAML to ensure strict adherence to the schema.

## EXECUTION PHASES

### [PHASE 1: Directory Mapping]
- Use the **Glob** tool to map the major folders and source files.
- Map them to `Module` dimension nodes in the hypergraph.

### [PHASE 2: Atomic Extraction]
- For each major file, use the **Read** tool to identify core functions or classes.
- Map these to `Atomic` dimension nodes.
- Link them to their parent Module via the `implements` edge.

### [PHASE 3: Edge Detection]
- Trace the data flow. If File A imports File B, File B is a dependency.
- Map this to the `depends_on` edge.

## OUTPUT

- **If creating from scratch:** Use the **Write** tool to save the compiled YAML to `spec/compiled/architecture.yml`. Set all newly discovered nodes to `status: clean`.
- **If appending/modifying an existing graph:** Use the **Edit** tool to update `spec/compiled/architecture.yml`, then use the **Bash** tool to propagate the blast radius:
  ```bash
  python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [modified_node_ids]
  ```

Output a summary of the system's topological density: number of System, Module, and Atomic nodes discovered.
