#!/usr/bin/env python3
"""
hypergraph_updater.py — Hypergraph blast-radius propagation and provenance staging.

Subcommands
-----------
Blast-radius (original behavior, default when no flag given):
    python hypergraph_updater.py <path_to_yaml> <dirty_node_1> [dirty_node_2 ...]

Provenance staging:
    python hypergraph_updater.py --write-provenance <branch_name> \
        --superprd <name> --miniprd <path> --api-call-id <id>

    python hypergraph_updater.py --merge-provenance <branch_name> <yaml_path>

    python hypergraph_updater.py --cleanup-provenance <branch_name>
"""

import fcntl
import os
import re
import sys
from datetime import datetime, timezone
from typing import List, Optional, Set

import yaml


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROVENANCE_STAGING_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", ".provenance_staging"
)
PROVENANCE_STAGING_DIR = os.path.normpath(PROVENANCE_STAGING_DIR)


# ---------------------------------------------------------------------------
# Original function — preserved exactly
# ---------------------------------------------------------------------------

def propagate_blast_radius(yaml_path: str, dirty_node_ids: List[str]) -> None:
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

    # Task 9 (Provenance MiniPRD) — exclusive lock on all writes to architecture.yml
    with open(yaml_path, 'w') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yaml.dump(data, f, sort_keys=False, default_flow_style=False)
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)

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


# ---------------------------------------------------------------------------
# Task 4 & 5 — Provenance staging write (with exclusive file lock)
# ---------------------------------------------------------------------------

def write_provenance_staging(
    branch_name: str,
    superprd_name: str,
    miniprd_path: str,
    api_call_id: str,
) -> None:
    """
    Write provenance staging entry to .provenance_staging/<branch_name>.yml.

    Called after daemon exit 0 (success) and BEFORE any rebase attempt.
    Acquires an exclusive lock before writing to prevent TOCTOU corruption.
    """
    os.makedirs(PROVENANCE_STAGING_DIR, exist_ok=True)

    staging_file = os.path.join(PROVENANCE_STAGING_DIR, f"{branch_name}.yml")
    iso_now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    payload = {
        "superprd_name": superprd_name,
        "branch_name": branch_name,
        "timestamp": iso_now,
        "miniprd_path": miniprd_path,
        "api_call_id": api_call_id,
    }

    # Task 5 — exclusive lock before write
    with open(staging_file, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yaml.dump(payload, fh, sort_keys=False, default_flow_style=False)
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)

    print("=" * 50)
    print("PROVENANCE STAGING: WRITTEN")
    print("=" * 50)
    print(f"File:         {staging_file}")
    print(f"Branch:       {branch_name}")
    print(f"SuperPRD:     {superprd_name}")
    print(f"MiniPRD path: {miniprd_path}")
    print(f"API call ID:  {api_call_id}")
    print(f"Timestamp:    {iso_now}")
    print("=" * 50)


# ---------------------------------------------------------------------------
# Task 6 — Rebase pipeline provenance merge
# ---------------------------------------------------------------------------

def _branch_name_to_node_id(branch_name: str) -> Optional[str]:
    """
    Derive the architecture node ID from a branch name.

    Convention: branch name contains the node ID with hyphens in place of
    underscores.  E.g. 'feature/hyper-orchestrator' → 'hyper_orchestrator'.

    Strategy:
      1. Strip common prefixes (feature/, bugfix/, fix/, chore/, hotfix/).
      2. Replace hyphens with underscores.
      3. Return the result as the candidate node ID.
    """
    # Strip common branch prefixes
    stripped = re.sub(r"^(feature|bugfix|fix|chore|hotfix)/", "", branch_name)
    # Replace hyphens with underscores to match node ID convention
    node_id = stripped.replace("-", "_")
    return node_id if node_id else None


def merge_provenance_to_arch(branch_name: str, yaml_path: str) -> None:
    """
    Read .provenance_staging/<branch_name>.yml and write the _provenance block
    to the relevant node in architecture.yml.

    The relevant node is identified by converting the branch name to a node ID
    using _branch_name_to_node_id().  Acquires an exclusive lock on
    architecture.yml during write.
    """
    staging_file = os.path.join(PROVENANCE_STAGING_DIR, f"{branch_name}.yml")

    if not os.path.exists(staging_file):
        print(
            f"ERROR: Provenance staging file not found: {staging_file}\n"
            f"Run --write-provenance first."
        )
        sys.exit(1)

    # Read staging file (shared lock)
    with open(staging_file, "r") as fh:
        fcntl.flock(fh, fcntl.LOCK_SH)
        try:
            staging_data = yaml.safe_load(fh) or {}
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)

    node_id = _branch_name_to_node_id(branch_name)
    if not node_id:
        print(f"ERROR: Could not derive node ID from branch name: '{branch_name}'")
        sys.exit(1)

    if not os.path.exists(yaml_path):
        print(f"ERROR: architecture YAML not found: {yaml_path}")
        sys.exit(1)

    # Read architecture.yml (shared lock)
    with open(yaml_path, "r") as fh:
        fcntl.flock(fh, fcntl.LOCK_SH)
        try:
            arch_data = yaml.safe_load(fh) or {"nodes": []}
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)

    nodes = arch_data.get("nodes", []) or []
    matched = False

    for node in nodes:
        if node.get("id") == node_id:
            node["_provenance"] = {
                "superprd_name": staging_data.get("superprd_name", ""),
                "branch_name": staging_data.get("branch_name", branch_name),
                "timestamp": staging_data.get("timestamp", ""),
                "miniprd_path": staging_data.get("miniprd_path", ""),
                "api_call_id": staging_data.get("api_call_id", ""),
            }
            matched = True
            break

    if not matched:
        print(
            f"WARNING: Node ID '{node_id}' (derived from branch '{branch_name}') "
            f"not found in {yaml_path}. Provenance block not written."
        )
        # Exit non-fatally — the node may not yet exist during early feature work.
        sys.exit(0)

    # Write back with exclusive lock (Task 6)
    with open(yaml_path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yaml.dump(arch_data, fh, sort_keys=False, default_flow_style=False)
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)

    print("=" * 50)
    print("PROVENANCE MERGE: SUCCESS")
    print("=" * 50)
    print(f"Branch:    {branch_name}")
    print(f"Node ID:   {node_id}")
    print(f"YAML:      {yaml_path}")
    print("=" * 50)


# ---------------------------------------------------------------------------
# Task 7 — Staging file cleanup
# ---------------------------------------------------------------------------

def cleanup_provenance_staging(branch_name: str) -> None:
    """Delete .provenance_staging/<branch_name>.yml after provenance is merged."""
    staging_file = os.path.join(PROVENANCE_STAGING_DIR, f"{branch_name}.yml")

    if os.path.exists(staging_file):
        os.remove(staging_file)
        print(f"[cleanup] Removed provenance staging file: {staging_file}")
    else:
        print(f"[cleanup] Staging file not found (already removed?): {staging_file}")

    # Attempt directory cleanup
    cleanup_staging_dir()


def cleanup_staging_dir() -> None:
    """Remove .provenance_staging/ directory if it is empty (ignoring .gitkeep)."""
    if not os.path.isdir(PROVENANCE_STAGING_DIR):
        return

    contents = os.listdir(PROVENANCE_STAGING_DIR)
    # Ignore .gitkeep placeholder
    meaningful = [f for f in contents if f != ".gitkeep"]

    if not meaningful:
        # Only .gitkeep (or truly empty) — leave the directory in place so
        # git tracks the placeholder.  If there is no .gitkeep either, we can
        # safely remove the directory.
        if not contents:
            os.rmdir(PROVENANCE_STAGING_DIR)
            print(f"[cleanup] Removed empty directory: {PROVENANCE_STAGING_DIR}")
        else:
            print(
                f"[cleanup] Directory contains only placeholder files; leaving in place: "
                f"{PROVENANCE_STAGING_DIR}"
            )
    else:
        print(
            f"[cleanup] Directory not empty ({len(meaningful)} file(s) remain); "
            f"skipping removal: {PROVENANCE_STAGING_DIR}"
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_and_dispatch(argv: List[str]) -> None:
    """Parse CLI arguments and dispatch to the appropriate function."""

    if not argv:
        _print_usage()
        sys.exit(1)

    first = argv[0]

    # ------------------------------------------------------------------
    # --write-provenance <branch_name> --superprd <name>
    #     --miniprd <path> --api-call-id <id>
    # ------------------------------------------------------------------
    if first == "--write-provenance":
        if len(argv) < 8:
            print(
                "Usage: hypergraph_updater.py --write-provenance <branch_name> "
                "--superprd <name> --miniprd <path> --api-call-id <id>"
            )
            sys.exit(1)

        branch_name = argv[1]
        args = argv[2:]

        def _flag_value(flag: str) -> Optional[str]:
            try:
                idx = args.index(flag)
                return args[idx + 1]
            except (ValueError, IndexError):
                return None

        superprd = _flag_value("--superprd")
        miniprd = _flag_value("--miniprd")
        api_call_id = _flag_value("--api-call-id")

        missing = [f for f, v in [("--superprd", superprd), ("--miniprd", miniprd), ("--api-call-id", api_call_id)] if v is None]
        if missing:
            print(f"ERROR: Missing required flag(s): {', '.join(missing)}")
            sys.exit(1)

        write_provenance_staging(branch_name, superprd, miniprd, api_call_id)

    # ------------------------------------------------------------------
    # --merge-provenance <branch_name> <yaml_path>
    # ------------------------------------------------------------------
    elif first == "--merge-provenance":
        if len(argv) < 3:
            print(
                "Usage: hypergraph_updater.py --merge-provenance <branch_name> <yaml_path>"
            )
            sys.exit(1)
        merge_provenance_to_arch(argv[1], argv[2])

    # ------------------------------------------------------------------
    # --cleanup-provenance <branch_name>
    # ------------------------------------------------------------------
    elif first == "--cleanup-provenance":
        if len(argv) < 2:
            print(
                "Usage: hypergraph_updater.py --cleanup-provenance <branch_name>"
            )
            sys.exit(1)
        cleanup_provenance_staging(argv[1])

    # ------------------------------------------------------------------
    # Legacy blast-radius: <yaml_path> <dirty_node_1> [dirty_node_2 ...]
    # ------------------------------------------------------------------
    else:
        if len(argv) < 2:
            _print_usage()
            sys.exit(1)
        propagate_blast_radius(argv[0], argv[1:])


def _print_usage() -> None:
    print(
        "Usage:\n"
        "  Blast-radius propagation (original):\n"
        "    python hypergraph_updater.py <path_to_yaml> <dirty_node_1> [dirty_node_2 ...]\n\n"
        "  Provenance staging:\n"
        "    python hypergraph_updater.py --write-provenance <branch_name> "
        "--superprd <name> --miniprd <path> --api-call-id <id>\n"
        "    python hypergraph_updater.py --merge-provenance <branch_name> <yaml_path>\n"
        "    python hypergraph_updater.py --cleanup-provenance <branch_name>\n"
    )


if __name__ == "__main__":
    _parse_and_dispatch(sys.argv[1:])
