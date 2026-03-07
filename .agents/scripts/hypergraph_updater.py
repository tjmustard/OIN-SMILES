#!/usr/bin/env python3
import yaml
import sys
import os
from typing import List, Set

def propagate_blast_radius(yaml_path: str, dirty_node_ids: List[str]):
    try:
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f) or {'nodes': []}
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Hypergraph not found at {yaml_path}")
        sys.exit(1)

    nodes = {node['id']: node for node in data.get('nodes', [])}
    
    for n_id in dirty_node_ids:
        if n_id not in nodes:
            print(f"FATAL: Dirty node '{n_id}' does not exist in hypergraph.")
            sys.exit(1)
            
    for n_id in dirty_node_ids:
        nodes[n_id]['status'] = 'dirty'
        
    queue = list(dirty_node_ids)
    processed: Set[str] = set(dirty_node_ids)

    while queue:
        current_id = queue.pop(0)
        current_node = nodes[current_id]
        edges = current_node.get('edges', {})

        for parent_id in edges.get('implements', []):
            if parent_id in nodes and parent_id not in processed:
                nodes[parent_id]['status'] = 'needs_review'
                processed.add(parent_id)
                queue.append(parent_id)
                
        for node_id, node_data in nodes.items():
            if current_id in node_data.get('edges', {}).get('depends_on', []):
                if node_id not in processed:
                    nodes[node_id]['status'] = 'needs_review'
                    processed.add(node_id)
                    queue.append(node_id)

    with open(yaml_path, 'w') as f:
        yaml.dump(data, f, sort_keys=False, default_flow_style=False)
        
    print("=" * 50)
    print("HYPERGRAPH UPDATER: SUCCESS")
    print("=" * 50)
    print(f"Target YAML: {yaml_path}")
    print(f"Nodes dirtied by direct execution: {dirty_node_ids}")
    print(f"Total nodes needing review (Blast Radius): {list(processed)}")
    print("=" * 50)
    print("AGENT INSTRUCTION: If you are the Builder Agent, your job is complete.")
    print("If you are the Auditor Agent, you must now semantically rewrite the")
    print("'outputs', 'inputs', and 'description' of these nodes based on the")
    print("new code, and return their state to 'clean'.")
    print("=" * 50)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python hypergraph_updater.py <path_to_yaml> <dirty_node_1> [dirty_node_2 ...]")
        sys.exit(1)
    propagate_blast_radius(sys.argv[1], sys.argv[2:])