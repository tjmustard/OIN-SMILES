#!/usr/bin/env python3
"""
semantic_graph_merger.py — Deterministic semantic merger for architecture.yml conflicts.

Called during git rebase when architecture.yml (or architecture.yaml) has conflicts.
NEVER uses AI. All merge decisions are deterministic.

CLI:
    python semantic_graph_merger.py <main_yaml> <branch_yaml> <output_yaml>
"""

import fcntl
import os
import re
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

import yaml


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FAILED_WORKFLOWS_PATH = os.path.normpath(os.path.join(
    os.path.dirname(__file__), "..", "..", "spec", "active", "FAILED_WORKFLOWS.md"
))
FAILED_WORKFLOWS_ARCHIVE_DIR = os.path.normpath(os.path.join(
    os.path.dirname(__file__), "..", "..", "spec", "archive"
))

REQUIRED_NODE_FIELDS = {"id", "dimension", "status"}


# ---------------------------------------------------------------------------
# Task 8 — File-type routing guard
# ---------------------------------------------------------------------------

def _is_valid_architecture_file(filepath: str) -> bool:
    """
    Returns True if the file is a valid architecture YAML target.

    Criteria (OR):
      (a) filename (basename without extension) contains 'architecture' (case-insensitive)
      (b) file extension is '.yml' or '.yaml'
    """
    basename = os.path.basename(filepath)
    name_without_ext, ext = os.path.splitext(basename)

    name_matches = "architecture" in name_without_ext.lower()
    ext_matches = ext.lower() in (".yml", ".yaml")

    return name_matches or ext_matches


def _guard_file_type(filepath: str) -> None:
    """Exits with error if the file is not a valid architecture/YAML target."""
    if not _is_valid_architecture_file(filepath):
        print(
            f"ERROR: semantic_graph_merger only accepts architecture YAML files.\n"
            f"  Rejected: {filepath}\n"
            f"  File must contain 'architecture' in its name or have a .yml/.yaml extension."
        )
        sys.exit(1)


# ---------------------------------------------------------------------------
# Task 11 — Extension normalization (canonical path resolution)
# ---------------------------------------------------------------------------

def _canonical_yaml_path(filepath: str) -> str:
    """
    Normalize .yml / .yaml paths to a single canonical path.

    Preference order:
      1. Use the path exactly as given if it exists on disk.
      2. If the .yml variant exists on disk, return that.
      3. If the .yaml variant exists on disk, return that.
      4. Return the path as-is (caller handles missing-file error).
    """
    if os.path.exists(filepath):
        return filepath

    root, ext = os.path.splitext(filepath)
    alt_ext = ".yaml" if ext.lower() == ".yml" else ".yml"
    alt_path = root + alt_ext

    if os.path.exists(alt_path):
        return alt_path

    return filepath


# ---------------------------------------------------------------------------
# Task 9 — Locked file I/O helpers
# ---------------------------------------------------------------------------

def _read_yaml_locked(filepath: str) -> Dict[str, Any]:
    """Read and parse a YAML file under a shared (read) lock."""
    with open(filepath, "r") as fh:
        fcntl.flock(fh, fcntl.LOCK_SH)
        try:
            data = yaml.safe_load(fh) or {"nodes": []}
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)
    return data


def _write_yaml_locked(filepath: str, data: Dict[str, Any]) -> None:
    """Write YAML data to a file under an exclusive lock."""
    with open(filepath, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yaml.dump(data, fh, sort_keys=False, default_flow_style=False)
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


# ---------------------------------------------------------------------------
# Task 13 — Pre-merge backup
# ---------------------------------------------------------------------------

def _backup_file(filepath: str) -> str:
    """Copy filepath to filepath.bak and return the backup path."""
    bak_path = filepath + ".bak"
    shutil.copy2(filepath, bak_path)
    return bak_path


def _restore_backup(filepath: str, bak_path: str) -> None:
    """Restore backup to original path."""
    if os.path.exists(bak_path):
        shutil.copy2(bak_path, filepath)


def _remove_backup(bak_path: str) -> None:
    """Remove backup file after successful validation."""
    if os.path.exists(bak_path):
        os.remove(bak_path)


# ---------------------------------------------------------------------------
# FAILED_WORKFLOWS.md cap + rotation helpers (Task 2 Provenance MiniPRD)
# ---------------------------------------------------------------------------

def _count_failed_headers(path: str) -> int:
    if not os.path.exists(path):
        return 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read().count("## [FAILED]")
    except OSError:
        return 0


def _rotate_failed_workflows() -> None:
    """Archive FAILED_WORKFLOWS.md when entry count reaches 50."""
    os.makedirs(FAILED_WORKFLOWS_ARCHIVE_DIR, exist_ok=True)
    iso = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = os.path.join(FAILED_WORKFLOWS_ARCHIVE_DIR, f"FAILED_WORKFLOWS_{iso}.md")
    shutil.copy2(FAILED_WORKFLOWS_PATH, dest)
    with open(FAILED_WORKFLOWS_PATH, "w", encoding="utf-8") as fh:
        fh.write("# Failed Workflows\n\n")
    print(f"[ROTATE] FAILED_WORKFLOWS.md archived to {dest}")


# ---------------------------------------------------------------------------
# FAILED_WORKFLOWS.md writer
# ---------------------------------------------------------------------------

def _append_failed_workflow(entry: str) -> None:
    """Append a failure entry to FAILED_WORKFLOWS.md with 50-entry cap + rotation."""
    os.makedirs(os.path.dirname(FAILED_WORKFLOWS_PATH), exist_ok=True)

    if not os.path.exists(FAILED_WORKFLOWS_PATH):
        header = (
            "# Failed Workflows\n\n"
            "This file is auto-maintained by `hyper_orchestrator.py`. "
            "Entries are appended on workflow failure.\n"
            "Maximum 50 entries — overflow is archived to "
            "`spec/archive/FAILED_WORKFLOWS_<timestamp>.md`.\n\n---\n"
        )
        with open(FAILED_WORKFLOWS_PATH, "w") as fh:
            fh.write(header)

    if _count_failed_headers(FAILED_WORKFLOWS_PATH) >= 50:
        _rotate_failed_workflows()

    with open(FAILED_WORKFLOWS_PATH, "a") as fh:
        fh.write("\n" + entry + "\n")


def _build_rebase_conflict_entry(
    yaml_filename: str,
    node_id: str,
    node_a: Dict[str, Any],
    node_b: Dict[str, Any],
) -> str:
    """Build a REBASE_CONFLICT markdown entry for FAILED_WORKFLOWS.md."""
    iso_now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    branch_context = os.environ.get("CURRENT_BRANCH", "unknown-branch")

    node_a_yaml = yaml.dump(node_a, sort_keys=False, default_flow_style=False).rstrip()
    node_b_yaml = yaml.dump(node_b, sort_keys=False, default_flow_style=False).rstrip()

    return (
        f"## [FAILED] {branch_context} — {iso_now}\n\n"
        f"- **Failure Type:** REBASE_CONFLICT\n"
        f"- **File:** {yaml_filename}\n"
        f"- **Conflicting Node ID:** {node_id}\n\n"
        f"### Branch A version:\n"
        f"```yaml\n{node_a_yaml}\n```\n\n"
        f"### Branch B version:\n"
        f"```yaml\n{node_b_yaml}\n```\n\n"
        f"---"
    )


def _build_schema_validation_failure_entry(
    yaml_filename: str,
    reason: str,
) -> str:
    """Build an ARCH_SCHEMA_VALIDATION_FAILURE markdown entry."""
    iso_now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    branch_context = os.environ.get("CURRENT_BRANCH", "unknown-branch")

    return (
        f"## [FAILED] {branch_context} — {iso_now}\n\n"
        f"- **Failure Type:** ARCH_SCHEMA_VALIDATION_FAILURE\n"
        f"- **File:** {yaml_filename}\n"
        f"- **Reason:** {reason}\n\n"
        f"---"
    )


# ---------------------------------------------------------------------------
# Task 10 — Merge algorithm with duplicate-ID conflict policy
# ---------------------------------------------------------------------------

def _nodes_are_identical(node_a: Dict[str, Any], node_b: Dict[str, Any]) -> bool:
    """Deep equality check for two node dicts."""
    return node_a == node_b


def _merge_node_lists(
    main_nodes: List[Dict[str, Any]],
    branch_nodes: List[Dict[str, Any]],
    yaml_filename: str,
) -> Tuple[Optional[List[Dict[str, Any]]], Optional[str]]:
    """
    Merge two node lists deterministically.

    Returns:
        (merged_list, None)          on success
        (None, conflict_entry_str)   when a REBASE_CONFLICT is detected
    """
    main_by_id: Dict[str, Dict[str, Any]] = {}
    for node in main_nodes:
        main_by_id[node["id"]] = node

    branch_by_id: Dict[str, Dict[str, Any]] = {}
    for node in branch_nodes:
        branch_by_id[node["id"]] = node

    merged_order: List[str] = []
    # Main-branch nodes first (preserve order)
    for node in main_nodes:
        if node["id"] not in merged_order:
            merged_order.append(node["id"])
    # New-only branch nodes appended
    for node in branch_nodes:
        if node["id"] not in merged_order:
            merged_order.append(node["id"])

    merged_by_id: Dict[str, Dict[str, Any]] = {}

    for node_id in merged_order:
        in_main = node_id in main_by_id
        in_branch = node_id in branch_by_id

        if in_main and not in_branch:
            # Present only in main — include as-is
            merged_by_id[node_id] = main_by_id[node_id]

        elif in_branch and not in_main:
            # Present only in branch — include as-is
            merged_by_id[node_id] = branch_by_id[node_id]

        else:
            # Present in both — apply Task 10 policy
            node_a = main_by_id[node_id]
            node_b = branch_by_id[node_id]

            if _nodes_are_identical(node_a, node_b):
                # Identical content → de-duplicate silently, keep one copy
                merged_by_id[node_id] = node_a
            else:
                # Differing content → abort; build conflict entry
                conflict_entry = _build_rebase_conflict_entry(
                    yaml_filename, node_id, node_a, node_b
                )
                return None, conflict_entry

    merged_list = [merged_by_id[nid] for nid in merged_order]
    return merged_list, None


# ---------------------------------------------------------------------------
# Task 12 — Post-merge schema validation
# ---------------------------------------------------------------------------

def _validate_architecture_yaml(filepath: str) -> Tuple[bool, str]:
    """
    Run post-merge schema validation checks.

    Returns (True, "") on success, or (False, reason_string) on failure.
    """
    # Check 1: parseable YAML
    try:
        with open(filepath, "r") as fh:
            fcntl.flock(fh, fcntl.LOCK_SH)
            try:
                data = yaml.safe_load(fh)
            finally:
                fcntl.flock(fh, fcntl.LOCK_UN)
    except yaml.YAMLError as exc:
        return False, f"YAML parse error: {exc}"

    if data is None:
        return False, "File parsed as empty/null YAML"

    nodes = data.get("nodes", [])

    # Check 2: no duplicate IDs
    seen_ids: List[str] = []
    for node in nodes:
        nid = node.get("id")
        if nid in seen_ids:
            return False, f"Duplicate node ID: '{nid}'"
        seen_ids.append(nid)

    all_ids = set(seen_ids)

    # Check 3: all depends_on references exist
    for node in nodes:
        edges = node.get("edges", {}) or {}
        deps = edges.get("depends_on", []) or []
        for dep_id in deps:
            if dep_id not in all_ids:
                return False, (
                    f"Node '{node.get('id')}' has depends_on reference to "
                    f"unknown node ID '{dep_id}'"
                )

    # Check 4: required fields present on every node
    for node in nodes:
        missing = REQUIRED_NODE_FIELDS - set(node.keys())
        if missing:
            return False, (
                f"Node '{node.get('id', '<unknown>')}' is missing required "
                f"field(s): {sorted(missing)}"
            )

    return True, ""


# ---------------------------------------------------------------------------
# Top-level merge entry point
# ---------------------------------------------------------------------------

def semantic_merge(main_yaml_raw: str, branch_yaml_raw: str, output_yaml_raw: str) -> None:
    """
    Perform deterministic semantic merge of two architecture YAML files.

    Steps:
      1. Guard file types (Task 8)
      2. Normalize extensions (Task 11)
      3. Take pre-merge backup (Task 13)
      4. Read both files under shared locks (Task 9)
      5. Merge node lists (Task 10)
      6. Write merged result under exclusive lock (Task 9)
      7. Post-merge schema validation (Task 12)
      8. Remove backup on success / restore on failure (Task 13)
    """
    # Task 8 — guard all three paths
    for path in (main_yaml_raw, branch_yaml_raw, output_yaml_raw):
        _guard_file_type(path)

    # Task 11 — normalize extensions to canonical on-disk paths
    main_yaml = _canonical_yaml_path(main_yaml_raw)
    branch_yaml = _canonical_yaml_path(branch_yaml_raw)
    output_yaml = _canonical_yaml_path(output_yaml_raw)

    yaml_filename = os.path.basename(output_yaml)

    # Validate input files exist
    for label, path in [("main_yaml", main_yaml), ("branch_yaml", branch_yaml)]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # Task 13 — pre-merge backup (of the output/main file if it exists)
    bak_path: Optional[str] = None
    if os.path.exists(output_yaml):
        bak_path = _backup_file(output_yaml)
        print(f"[backup] {output_yaml} → {bak_path}")

    # Read both YAML files under shared locks
    print(f"[read]   {main_yaml}")
    main_data = _read_yaml_locked(main_yaml)
    print(f"[read]   {branch_yaml}")
    branch_data = _read_yaml_locked(branch_yaml)

    main_nodes: List[Dict[str, Any]] = main_data.get("nodes", []) or []
    branch_nodes: List[Dict[str, Any]] = branch_data.get("nodes", []) or []

    # Task 10 — merge with conflict policy
    merged_nodes, conflict_entry = _merge_node_lists(
        main_nodes, branch_nodes, yaml_filename
    )

    if conflict_entry is not None:
        # REBASE_CONFLICT — abort, do NOT touch output_yaml
        _append_failed_workflow(conflict_entry)
        print(
            f"ERROR: REBASE_CONFLICT detected. "
            f"Merge aborted. architecture.yml left unmodified.\n"
            f"Conflict details written to: {FAILED_WORKFLOWS_PATH}"
        )
        # Leave backup in place so caller can inspect original state
        sys.exit(1)

    # Build merged data structure (preserve non-node top-level keys from main)
    merged_data: Dict[str, Any] = dict(main_data)
    merged_data["nodes"] = merged_nodes

    # Task 9 — write under exclusive lock
    print(f"[write]  {output_yaml}")
    _write_yaml_locked(output_yaml, merged_data)

    # Task 12 — post-merge schema validation
    valid, reason = _validate_architecture_yaml(output_yaml)
    if not valid:
        # Restore backup
        if bak_path:
            _restore_backup(output_yaml, bak_path)
            print(f"[restore] Backup restored: {bak_path} → {output_yaml}")

        failure_entry = _build_schema_validation_failure_entry(yaml_filename, reason)
        _append_failed_workflow(failure_entry)
        print(
            f"ERROR: ARCH_SCHEMA_VALIDATION_FAILURE: {reason}\n"
            f"Merged file rejected; original restored.\n"
            f"Failure details written to: {FAILED_WORKFLOWS_PATH}"
        )
        sys.exit(1)

    # Task 13 — validation passed, remove backup
    if bak_path:
        _remove_backup(bak_path)
        print(f"[cleanup] Backup removed: {bak_path}")

    node_count = len(merged_nodes)
    print("=" * 50)
    print("SEMANTIC GRAPH MERGER: SUCCESS")
    print("=" * 50)
    print(f"Main YAML:    {main_yaml}")
    print(f"Branch YAML:  {branch_yaml}")
    print(f"Output YAML:  {output_yaml}")
    print(f"Merged nodes: {node_count}")
    print("=" * 50)


# ---------------------------------------------------------------------------
# Git conflict-marker support (1-arg invocation from hyper_orchestrator.py)
# ---------------------------------------------------------------------------

def _parse_conflict_markers(content: str) -> tuple:
    """
    Split a git conflict-marker file into (ours, theirs) strings.

    Handles multiple conflict regions within the same file; non-conflicted
    lines are included in both output strings verbatim.
    """
    ours_lines: List[str] = []
    theirs_lines: List[str] = []
    in_conflict = False
    in_theirs_section = False

    for line in content.splitlines(keepends=True):
        stripped = line.rstrip("\n").rstrip("\r")
        if stripped.startswith("<<<<<<<"):
            in_conflict = True
            in_theirs_section = False
        elif stripped == "=======" and in_conflict:
            in_theirs_section = True
        elif stripped.startswith(">>>>>>>") and in_conflict:
            in_conflict = False
            in_theirs_section = False
        elif in_conflict and not in_theirs_section:
            ours_lines.append(line)
        elif in_conflict and in_theirs_section:
            theirs_lines.append(line)
        else:
            ours_lines.append(line)
            theirs_lines.append(line)

    return "".join(ours_lines), "".join(theirs_lines)


def _merge_from_conflict_markers(filepath: str) -> None:
    """
    1-arg mode: read a git-conflict-marked file, extract ours/theirs YAML
    versions, and call semantic_merge() with temporary intermediate files.
    """
    _guard_file_type(filepath)

    try:
        with open(filepath, "r", encoding="utf-8") as fh:
            content = fh.read()
    except OSError as exc:
        print(f"ERROR: Cannot read {filepath}: {exc}")
        sys.exit(1)

    if "<<<<<<<" not in content:
        print(f"ERROR: {filepath} contains no git conflict markers (<<<<<<<). "
              "Nothing to merge.")
        sys.exit(1)

    ours, theirs = _parse_conflict_markers(content)

    main_tmp = branch_tmp = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yml", delete=False,
            prefix="sgm_main_", encoding="utf-8"
        ) as f:
            f.write(ours)
            main_tmp = f.name
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yml", delete=False,
            prefix="sgm_branch_", encoding="utf-8"
        ) as f:
            f.write(theirs)
            branch_tmp = f.name

        semantic_merge(main_tmp, branch_tmp, filepath)
    finally:
        for tmp in (main_tmp, branch_tmp):
            if tmp:
                try:
                    os.unlink(tmp)
                except OSError:
                    pass


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) == 2:
        # 1-arg form: parse git conflict markers from the single file
        _merge_from_conflict_markers(sys.argv[1])
    elif len(sys.argv) == 4:
        # 3-arg form: explicit main / branch / output paths
        semantic_merge(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print(
            "Usage:\n"
            "  python semantic_graph_merger.py <conflicted_file>\n"
            "  python semantic_graph_merger.py <main_yaml> <branch_yaml> <output_yaml>"
        )
        sys.exit(1)
