"""
Integration tests for MiniPRD: Provenance Integration

Covers:
  - semantic_graph_merger.py routing guard (prefix + .yaml extension)
  - Duplicate node ID conflict policies (identical + differing)
  - Post-merge schema validation (duplicate IDs, missing depends_on)
  - Provenance staging lifecycle
  - FAILED_WORKFLOWS.md 50-entry cap via semantic_graph_merger
  - Concurrent write protection
  - Pre-merge backup restore on validation failure
"""

import fcntl
import os
import sys
import shutil
import threading
import time
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

# Ensure the scripts directory is importable
SCRIPTS_DIR = Path(__file__).parent.parent.parent / ".agents" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import semantic_graph_merger as sgm  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_yaml(nodes):
    """Return a YAML string for a simple architecture doc."""
    return yaml.dump({"nodes": nodes}, sort_keys=False, default_flow_style=False)


def _write(path, content):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(content, encoding="utf-8")
    return path


def _node(id_, description="desc", depends_on=None):
    n = {"id": id_, "dimension": "Module", "status": "clean", "description": description}
    n["edges"] = {"depends_on": depends_on or [], "implements": []}
    return n


# ---------------------------------------------------------------------------
# Test 1 — Routing guard: architecture prefix routes to merger
# ---------------------------------------------------------------------------

def test_routing_guard_architecture_prefix():
    """File named architecture_v2.yml should pass _is_valid_architecture_file."""
    assert sgm._is_valid_architecture_file("architecture_v2.yml") is True
    assert sgm._is_valid_architecture_file("ARCHITECTURE_GRAPH.YML") is True  # case-insensitive
    assert sgm._is_valid_architecture_file("config.py") is False               # negative case


# ---------------------------------------------------------------------------
# Test 2 — Routing guard: .yaml extension routes to merger
# ---------------------------------------------------------------------------

def test_routing_guard_yaml_extension():
    """Any .yaml file should pass the routing guard."""
    assert sgm._is_valid_architecture_file("architecture.yaml") is True
    assert sgm._is_valid_architecture_file("something.yaml") is True


# ---------------------------------------------------------------------------
# Test 3 — Duplicate node (identical): silent dedup
# ---------------------------------------------------------------------------

def test_duplicate_node_identical_dedup(tmp_path):
    node = _node("my_module", description="same")
    main_yaml = tmp_path / "main.yml"
    branch_yaml = tmp_path / "branch.yml"
    output_yaml = tmp_path / "architecture.yml"

    _write(main_yaml, _make_yaml([node]))
    _write(branch_yaml, _make_yaml([node]))
    _write(output_yaml, _make_yaml([node]))

    sgm.semantic_merge(str(main_yaml), str(branch_yaml), str(output_yaml))

    result = yaml.safe_load(output_yaml.read_text())
    assert len(result["nodes"]) == 1
    assert result["nodes"][0]["id"] == "my_module"


# ---------------------------------------------------------------------------
# Test 4 — Duplicate node (differing): REBASE_CONFLICT + output unmodified
# ---------------------------------------------------------------------------

def test_duplicate_node_differing_rebase_conflict(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Create FAILED_WORKFLOWS dir
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"

    node_a = _node("clash", description="version A")
    node_b = _node("clash", description="version B")
    main_yaml = tmp_path / "main.yml"
    branch_yaml = tmp_path / "branch.yml"
    output_yaml = tmp_path / "architecture.yml"
    original_content = _make_yaml([node_a])

    _write(main_yaml, _make_yaml([node_a]))
    _write(branch_yaml, _make_yaml([node_b]))
    _write(output_yaml, original_content)

    # Patch the FAILED_WORKFLOWS path inside the module
    monkeypatch.setattr(
        sgm,
        "FAILED_WORKFLOWS_PATH",
        str(fw_path),
    )

    with pytest.raises(SystemExit) as exc_info:
        sgm.semantic_merge(str(main_yaml), str(branch_yaml), str(output_yaml))

    assert exc_info.value.code == 1
    # Output file unmodified
    assert output_yaml.read_text(encoding="utf-8") == original_content
    # REBASE_CONFLICT entry written
    assert fw_path.exists()
    assert "REBASE_CONFLICT" in fw_path.read_text()


# ---------------------------------------------------------------------------
# Test 5 — Post-merge schema validation: duplicate IDs
# ---------------------------------------------------------------------------

def test_schema_validation_duplicate_ids(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(sgm, "FAILED_WORKFLOWS_PATH", str(fw_path))

    node_a = _node("alpha")
    node_b = _node("beta")
    node_dup = _node("alpha")  # duplicate id

    # The merger itself won't produce this; simulate by crafting branch_yaml
    # that when "merged" produces duplicate IDs. We do it by having both
    # main and branch each add a different node, then manually inject dup.
    # Simplest: mock _merge_node_lists to return a dup-id list.
    main_yaml = tmp_path / "main.yml"
    branch_yaml = tmp_path / "branch.yml"
    output_yaml = tmp_path / "architecture.yml"
    original = _make_yaml([node_a])

    _write(main_yaml, _make_yaml([node_a]))
    _write(branch_yaml, _make_yaml([node_b]))
    _write(output_yaml, original)

    # Patch _merge_node_lists to return a list with duplicate IDs
    def _fake_merge(main_nodes, branch_nodes, yaml_filename):
        return [node_a, node_dup], None  # 'alpha' appears twice

    monkeypatch.setattr(sgm, "_merge_node_lists", _fake_merge)

    with pytest.raises(SystemExit) as exc_info:
        sgm.semantic_merge(str(main_yaml), str(branch_yaml), str(output_yaml))

    assert exc_info.value.code == 1
    # Backup restored — output is byte-identical to pre-merge state
    assert output_yaml.read_text(encoding="utf-8") == original
    assert fw_path.exists()
    assert "ARCH_SCHEMA_VALIDATION_FAILURE" in fw_path.read_text()


# ---------------------------------------------------------------------------
# Test 6 — Post-merge schema validation: missing depends_on reference
# ---------------------------------------------------------------------------

def test_schema_validation_missing_depends_on(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(sgm, "FAILED_WORKFLOWS_PATH", str(fw_path))

    node_bad = _node("alpha", depends_on=["nonexistent_node"])
    main_yaml = tmp_path / "main.yml"
    branch_yaml = tmp_path / "branch.yml"
    output_yaml = tmp_path / "architecture.yml"
    original = _make_yaml([node_bad])

    _write(main_yaml, _make_yaml([node_bad]))
    _write(branch_yaml, _make_yaml([node_bad]))
    _write(output_yaml, original)

    with pytest.raises(SystemExit) as exc_info:
        sgm.semantic_merge(str(main_yaml), str(branch_yaml), str(output_yaml))

    assert exc_info.value.code == 1
    assert output_yaml.read_text(encoding="utf-8") == original
    assert fw_path.exists()
    assert "ARCH_SCHEMA_VALIDATION_FAILURE" in fw_path.read_text()


# ---------------------------------------------------------------------------
# Test 7 — Provenance staging lifecycle
# ---------------------------------------------------------------------------

def test_provenance_staging_lifecycle(tmp_path):
    """
    Verify that _write_provenance_staging() (from hyper_orchestrator) creates
    the staging file with the required payload fields, and that
    _merge_provenance_into_arch() moves them into architecture.yml and
    deletes the staging file.
    """
    # Import orchestrator helpers
    orch_scripts = SCRIPTS_DIR
    sys.path.insert(0, str(orch_scripts))
    import hyper_orchestrator as orch

    staging_dir = tmp_path / ".provenance_staging"
    arch_path = tmp_path / "architecture.yml"

    # Write a minimal architecture.yml with a target node
    node = {"id": "provenance_integration", "dimension": "Module", "status": "dirty"}
    arch_path.write_text(
        yaml.dump({"nodes": [node]}, sort_keys=False), encoding="utf-8"
    )

    # Patch constants in orchestrator module
    original_staging = orch.PROVENANCE_STAGING_DIR
    original_arch = orch.ARCH_YAML_PATH
    orch.PROVENANCE_STAGING_DIR = staging_dir
    orch.ARCH_YAML_PATH = arch_path

    try:
        branch = "feature/miniprd-provenance-integration"
        miniprd_path = tmp_path / "MiniPRD_Provenance_Integration.md"
        miniprd_path.write_text("# MiniPRD", encoding="utf-8")

        # Step 1: daemon success → staging file written
        orch._write_provenance_staging(
            branch_name=branch,
            superprd_name="SuperPRD_Test",
            miniprd_path=miniprd_path,
            message_id="msg_test_123",
        )
        safe_branch = branch.replace("/", "_")
        staging_file = staging_dir / f"{safe_branch}.yml"
        assert staging_file.exists(), "Staging file must exist after daemon success"

        payload = yaml.safe_load(staging_file.read_text())
        assert payload["superprd_name"] == "SuperPRD_Test"
        assert payload["branch_name"] == branch
        assert payload["api_call_id"] == "msg_test_123"
        assert "timestamp" in payload

        # Step 2: rebase success → _provenance written into arch, staging deleted
        orch._merge_provenance_into_arch(branch)

        data = yaml.safe_load(arch_path.read_text())
        prov_node = next(
            (n for n in data["nodes"] if n["id"] == "provenance_integration"), None
        )
        assert prov_node is not None, "Node must exist in architecture.yml"
        assert "_provenance" in prov_node, "_provenance block must be written"
        assert prov_node["_provenance"]["api_call_id"] == "msg_test_123"
        assert not staging_file.exists(), "Staging file must be deleted after merge"
    finally:
        orch.PROVENANCE_STAGING_DIR = original_staging
        orch.ARCH_YAML_PATH = original_arch


# ---------------------------------------------------------------------------
# Test 8 — FAILED_WORKFLOWS.md 50-entry cap via semantic_graph_merger
# ---------------------------------------------------------------------------

def test_failed_workflows_50_entry_cap(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    archive_dir = tmp_path / "spec" / "archive"
    fw_path.parent.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(sgm, "FAILED_WORKFLOWS_PATH", str(fw_path))
    monkeypatch.setattr(sgm, "FAILED_WORKFLOWS_ARCHIVE_DIR", str(archive_dir))

    # Populate with exactly 50 entries
    entries = "# Failed Workflows\n\n"
    for i in range(50):
        entries += f"\n## [FAILED] branch-{i} — 2025-01-01T00:00:00Z\n\n- entry {i}\n\n---\n"
    fw_path.write_text(entries, encoding="utf-8")

    # Trigger one more failure via semantic_graph_merger
    node_a = _node("x", description="v1")
    node_b = _node("x", description="v2")
    main_yaml = tmp_path / "main.yml"
    branch_yaml = tmp_path / "branch.yml"
    output_yaml = tmp_path / "architecture.yml"
    _write(main_yaml, _make_yaml([node_a]))
    _write(branch_yaml, _make_yaml([node_b]))
    _write(output_yaml, _make_yaml([node_a]))

    with pytest.raises(SystemExit):
        sgm.semantic_merge(str(main_yaml), str(branch_yaml), str(output_yaml))

    # Old file archived
    archived = list(archive_dir.glob("FAILED_WORKFLOWS_*.md"))
    assert len(archived) == 1, "Exactly one archive file must be created"
    assert "## [FAILED]" in archived[0].read_text()

    # Fresh file has exactly 1 entry
    assert fw_path.read_text().count("## [FAILED]") == 1


# ---------------------------------------------------------------------------
# Test 9 — Concurrent write protection
# ---------------------------------------------------------------------------

def test_concurrent_write_protection(tmp_path):
    """
    Two threads writing to architecture.yml via _write_yaml_locked() must
    not corrupt the file; the final result must be valid YAML.
    """
    arch = tmp_path / "architecture.yml"
    node = {"id": "a", "dimension": "Module", "status": "clean", "edges": {"depends_on": [], "implements": []}}
    arch.write_text(yaml.dump({"nodes": [node]}), encoding="utf-8")

    errors: list[Exception] = []

    def writer(node_id: str):
        try:
            data = {"nodes": [{"id": node_id, "dimension": "Module", "status": "clean",
                                "edges": {"depends_on": [], "implements": []}}]}
            sgm._write_yaml_locked(str(arch), data)
        except Exception as exc:
            errors.append(exc)

    t1 = threading.Thread(target=writer, args=("thread_a",))
    t2 = threading.Thread(target=writer, args=("thread_b",))
    t1.start()
    t2.start()
    t1.join(timeout=10)
    t2.join(timeout=10)

    assert not errors, f"Concurrent write raised: {errors}"
    # Final file must parse as valid YAML
    result = yaml.safe_load(arch.read_text())
    assert "nodes" in result
    assert len(result["nodes"]) >= 1


# ---------------------------------------------------------------------------
# Test 10 — Pre-merge backup restore on validation failure
# ---------------------------------------------------------------------------

def test_pre_merge_backup_restore(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(sgm, "FAILED_WORKFLOWS_PATH", str(fw_path))

    node_a = _node("stable")
    main_yaml = tmp_path / "main.yml"
    branch_yaml = tmp_path / "branch.yml"
    output_yaml = tmp_path / "architecture.yml"
    pre_merge_content = _make_yaml([node_a])

    _write(main_yaml, _make_yaml([node_a]))
    _write(branch_yaml, _make_yaml([node_a]))
    _write(output_yaml, pre_merge_content)

    # Force validation failure by injecting duplicate IDs after merge
    def _bad_merge(main_nodes, branch_nodes, yaml_filename):
        return [_node("stable"), _node("stable")], None  # duplicate

    monkeypatch.setattr(sgm, "_merge_node_lists", _bad_merge)

    with pytest.raises(SystemExit) as exc_info:
        sgm.semantic_merge(str(main_yaml), str(branch_yaml), str(output_yaml))

    assert exc_info.value.code == 1
    # architecture.yml byte-identical to pre-merge state (backup restored)
    assert output_yaml.read_text(encoding="utf-8") == pre_merge_content
    # MiniPRD Task 13: ".bak left in place on failure for forensic review"
    assert Path(str(output_yaml) + ".bak").exists()
