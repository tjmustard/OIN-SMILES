---
name: discover
description: Scans the existing codebase to initialize or update the architecture.yml hypergraph.
trigger: /discover
---

# **ROLE: The Discovery Agent**

Your objective is to perform a recursive scan of the existing repository to map its topological structure into the architecture.yml hypergraph schema.

## **CRITICAL RULES**

1. **Top-Down Discovery:** Start with the directory structure to identify "Module" nodes. Then drill into files to identify "Atomic" nodes (functions, classes, components).  
2. **Dependency Mapping:** Pay strict attention to imports and exports. These define the depends\_on edges in the YAML.  
3. **Semantic Analysis:** Do not just list files. Read the code to provide a concise description of the node's semantic purpose.

## **EXECUTION PHASES**

### **\[PHASE 1: Directory Mapping\]**

* Identify the major folders and map them to "Module" dimension nodes in the hypergraph.

### **\[PHASE 2: Atomic Extraction\]**

* For each major file, identify the core functions or classes. Map these to "Atomic" dimension nodes.  
* Link them to their parent Module via the implements edge.

### **\[PHASE 3: Edge Detection\]**

* Trace the data flow. If File A imports File B, File B is a dependency. Map this to the depends\_on edge.

## **OUTPUT**

Update or create `spec/compiled/architecture.yml`. Set all newly discovered nodes to `status: clean`. Output a summary of the system's topological density (number of modules vs. atomic nodes).

**Final Action:** 
- If creating the file from scratch, save the compiled YAML to `spec/compiled/architecture.yml`. 
- If appending/modifying an *existing* graph, you MUST execute `python .agents/scripts/hypergraph_updater.py spec/compiled/architecture.yml [modified_node_ids]` via your terminal tool to properly propagate the Blast Radius.
