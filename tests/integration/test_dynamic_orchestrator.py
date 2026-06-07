"""
Integration tests for MiniPRD: Dynamic Orchestrator

Covers hyper_orchestrator.py and hyper_daemon.py:
  - Global mutex prevents concurrent invocations
  - Startup validation rejects null-byte .md files
  - Branch collision detection exits before thread spawn
  - TOKEN_OVERRUN sentinel (EXIT_TOKEN_OVERRUN=3) breaks daemon loop
  - NO_TESTS_COLLECTED abort without incrementing iteration
  - Hung worker marked TIMEOUT_FAILURE after 600s wait
  - SIGTERM releases lock and writes partial FAILED_WORKFLOWS.md
  - Idempotency skips branches with existing _provenance
"""

import concurrent.futures
import fcntl
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch, call

import pytest
import yaml

SCRIPTS_DIR = Path(__file__).parent.parent.parent / ".agents" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import hyper_orchestrator as orch  # noqa: E402
import hyper_daemon as daemon       # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_fw_path(tmp_path):
    p = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _patch_fw(monkeypatch, tmp_path):
    fw = _make_fw_path(tmp_path)
    monkeypatch.setattr(orch, "FAILED_WORKFLOWS_PATH", fw)
    monkeypatch.setattr(daemon, "FAILED_WORKFLOWS_PATH", fw)
    return fw


# ---------------------------------------------------------------------------
# Test 1 — Global mutex: second process exits with clear error
# ---------------------------------------------------------------------------

def test_mutex_second_process_exits(tmp_path, monkeypatch):
    """
    If .orchestrator.lock is already held, orchestrate() must exit with an
    error before spawning any threads.
    """
    monkeypatch.chdir(tmp_path)
    lock_path = tmp_path / ".orchestrator.lock"
    fw = _patch_fw(monkeypatch, tmp_path)
    monkeypatch.setattr(orch, "LOCK_FILE", lock_path)
    monkeypatch.setattr(orch, "PID_LOCK_FILE", tmp_path / ".orchestrator.pid.lock")
    monkeypatch.setattr(orch, "LOG_DIR", tmp_path / ".agents" / "logs")

    # Hold the lock externally
    lock_fd = open(lock_path, "w")
    fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    try:
        t0 = time.monotonic()
        with pytest.raises(SystemExit) as exc_info:
            orch._acquire_mutex()
        elapsed = time.monotonic() - t0
        assert exc_info.value.code != 0
        assert "orchestrator" in str(exc_info.value).lower() or exc_info.value.code != 0
        assert elapsed < 2.0, f"Mutex rejection took {elapsed:.2f}s (must be < 2s)"
    finally:
        fcntl.flock(lock_fd, fcntl.LOCK_UN)
        lock_fd.close()


# ---------------------------------------------------------------------------
# Test 2 — Startup validation: null byte in .md file exits before threads
# ---------------------------------------------------------------------------

def test_startup_validation_null_byte(tmp_path, monkeypatch):
    """A .md file containing \\x00 must cause exit before any thread is spawned."""
    monkeypatch.chdir(tmp_path)
    superprd_dir = tmp_path / "specs"
    superprd_dir.mkdir()
    bad_md = superprd_dir / "bad_spec.md"
    # \x80\x81\x82 are invalid UTF-8 leading bytes — triggers UnicodeDecodeError
    bad_md.write_bytes(b"# Spec\n\x80\x81\x82")

    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")

    with pytest.raises(SystemExit) as exc_info:
        orch._validate_startup(superprd_dir, [bad_md], max_api_calls=50)

    assert exc_info.value.code != 0
    assert "UTF-8" in str(exc_info.value)


# ---------------------------------------------------------------------------
# Test 3 — Branch collision exits before spawning threads
# ---------------------------------------------------------------------------

def test_branch_collision_exits(tmp_path, monkeypatch):
    """If a target branch already exists, _validate_startup must exit listing it."""
    monkeypatch.chdir(tmp_path)
    superprd_dir = tmp_path / "specs"
    superprd_dir.mkdir()
    miniprd = superprd_dir / "MiniPRD_Dynamic_Orchestrator.md"
    miniprd.write_text("# MiniPRD\n", encoding="utf-8")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")

    expected_branch = orch._branch_name_for_miniprd(miniprd)

    def _fake_subprocess_run(cmd, **kwargs):
        result = SimpleNamespace()
        result.stdout = ""
        result.stderr = ""
        if isinstance(cmd, list) and "show-ref" in cmd:
            # Any branch ref lookup returns 0 (branch exists)
            result.returncode = 0
        elif isinstance(cmd, list) and cmd[:1] == ["du"]:
            result.returncode = 0
            result.stdout = "1000000\t."
        else:
            result.returncode = 1
        return result

    monkeypatch.setattr(orch.subprocess, "run", _fake_subprocess_run)

    with pytest.raises(SystemExit) as exc_info:
        orch._validate_startup(superprd_dir, [miniprd], max_api_calls=50)

    assert exc_info.value.code != 0
    # The exit message should contain the colliding branch name
    assert any(expected_branch in str(a) for a in exc_info.value.args)


# ---------------------------------------------------------------------------
# Test 4 — TOKEN_OVERRUN sentinel breaks daemon loop without counting iteration
# ---------------------------------------------------------------------------

def test_token_overrun_exits_loop(tmp_path, monkeypatch):
    """
    When hyper_fix.py exits with EXIT_TOKEN_OVERRUN (3), run_loop() must
    break immediately and return False without incrementing the iteration counter
    past the point of failure.
    """
    monkeypatch.chdir(tmp_path)
    fw = _patch_fw(monkeypatch, tmp_path)
    miniprd = tmp_path / "MiniPRD_Test.md"
    miniprd.write_text("# MiniPRD\n", encoding="utf-8")

    iteration_log = []

    def _fake_run_agent(script_name, args, env_vars=None):
        iteration_log.append(script_name)
        return True  # execute/audit succeed

    def _fake_pytest_collect():
        return 5  # tests found

    def _fake_oracle():
        return 1, "FAIL", ""  # oracle fails

    # Fix subprocess returns EXIT_TOKEN_OVERRUN
    def _fake_subprocess_run(cmd, **kwargs):
        r = SimpleNamespace()
        if "hyper_fix" in str(cmd):
            r.returncode = daemon.EXIT_TOKEN_OVERRUN
        else:
            r.returncode = 0
        return r

    monkeypatch.setattr(daemon, "_run_agent_script", _fake_run_agent)
    monkeypatch.setattr(daemon, "_pytest_collect_count", _fake_pytest_collect)
    monkeypatch.setattr(daemon, "_run_oracle_tests", _fake_oracle)
    monkeypatch.setattr(daemon.subprocess, "run", _fake_subprocess_run)
    monkeypatch.setattr(daemon, "_capture_baseline_interface", lambda p: {"foo"})

    result = daemon.run_loop(
        superprd_name="TestSuperPRD",
        miniprd_path=miniprd,
        orchestrated_mode=False,
    )

    assert result is False
    # Only one fix attempt was made (loop stopped after TOKEN_OVERRUN)
    # fix is called via subprocess, not _run_agent_script
    assert iteration_log.count("hyper_execute_runner") == 1
    # TOKEN_OVERRUN must be logged to FAILED_WORKFLOWS.md
    assert "TOKEN_OVERRUN" in fw.read_text(), "TOKEN_OVERRUN must be written to FAILED_WORKFLOWS.md"
    # Truncated output must not be written to disk; only the MiniPRD fixture file exists
    unexpected = [
        p for p in tmp_path.rglob("*")
        if p.is_file() and p.suffix in (".py", ".txt") and p != miniprd
    ]
    assert not unexpected, f"No code output should be written on TOKEN_OVERRUN; found: {unexpected}"


# ---------------------------------------------------------------------------
# Test 5 — NO_TESTS_COLLECTED does not increment iteration counter
# ---------------------------------------------------------------------------

def test_no_tests_collected_no_iteration(tmp_path, monkeypatch):
    """When pytest --collect-only returns 0, abort without counting it as an iteration."""
    monkeypatch.chdir(tmp_path)
    fw = _patch_fw(monkeypatch, tmp_path)
    miniprd = tmp_path / "MiniPRD_Test.md"
    miniprd.write_text("# MiniPRD\n", encoding="utf-8")

    monkeypatch.setattr(daemon, "_run_agent_script", lambda *a, **kw: True)
    monkeypatch.setattr(daemon, "_pytest_collect_count", lambda: 0)
    monkeypatch.setattr(daemon, "_capture_baseline_interface", lambda p: set())

    result = daemon.run_loop(
        superprd_name="Test",
        miniprd_path=miniprd,
        orchestrated_mode=False,
    )

    assert result is False
    assert fw.exists()
    fw_content = fw.read_text()
    assert "NO_TESTS_COLLECTED" in fw_content
    # iterations_attempted should be 0 (iteration - 1 = 1 - 1 = 0)
    assert "0/5" in fw_content


# ---------------------------------------------------------------------------
# Test 6 — Hung worker marked TIMEOUT_FAILURE after wait timeout
# ---------------------------------------------------------------------------

def test_hung_worker_timeout_failure(tmp_path, monkeypatch):
    """
    When concurrent.futures.wait times out, not_done futures must be marked
    TIMEOUT_FAILURE by the orchestrator — verified by driving orchestrate().
    """
    monkeypatch.chdir(tmp_path)
    fw = _patch_fw(monkeypatch, tmp_path)
    monkeypatch.setattr(orch, "LOCK_FILE", tmp_path / ".orchestrator.lock")
    monkeypatch.setattr(orch, "PID_LOCK_FILE", tmp_path / ".orchestrator.pid.lock")
    monkeypatch.setattr(orch, "LOG_DIR", tmp_path / ".agents" / "logs")
    monkeypatch.setattr(orch, "ARCH_YAML_PATH", tmp_path / "architecture.yml")
    monkeypatch.setattr(orch, "PROVENANCE_STAGING_DIR", tmp_path / ".provenance_staging")

    # Minimal file structure for pair loading (use "superprd" not "spec" —
    # _patch_fw already creates tmp_path/spec/active/)
    superprd_dir = tmp_path / "superprd"
    superprd_dir.mkdir()
    superprd_file = superprd_dir / "SuperPRD.md"
    superprd_file.write_text("# SuperPRD\n", encoding="utf-8")
    miniprd = superprd_dir / "MiniPRD_HungBranch.md"
    miniprd.write_text("# MiniPRD\n", encoding="utf-8")

    # A Future that is never resolved — simulates a hung worker
    hung_future = concurrent.futures.Future()

    # Mock out infrastructure to avoid real git / API / disk validation
    monkeypatch.setattr(orch, "_validate_startup", lambda *a, **k: None)
    monkeypatch.setattr(orch, "_write_pid_manifest", lambda *a, **k: None)
    monkeypatch.setattr(orch, "_remove_pid_manifest", lambda: None)
    monkeypatch.setattr(orch, "warn_stale_failed_branches", lambda: None)
    monkeypatch.setattr(
        orch, "_load_superprd_miniprd_pairs",
        lambda d: [(superprd_file, miniprd, "TestSuperPRD")],
    )

    # Return a real file so fcntl in orchestrate()'s finally block works
    lock_file = tmp_path / ".orchestrator.lock"
    lock_fd = open(lock_file, "w")
    monkeypatch.setattr(orch, "_acquire_mutex", lambda: lock_fd)

    # Executor returns the hung future; wait immediately puts it in not_done
    mock_executor = MagicMock()
    mock_executor.submit.return_value = hung_future
    mock_executor.shutdown = MagicMock()

    monkeypatch.setattr(
        concurrent.futures, "wait",
        lambda fs, timeout=None: (frozenset(), frozenset(fs)),
    )

    with patch("concurrent.futures.ThreadPoolExecutor", return_value=mock_executor):
        orch.orchestrate(superprd_dir=superprd_dir)

    # orchestrate()'s not_done loop must have written TIMEOUT_FAILURE
    assert fw.exists()
    fw_content = fw.read_text()
    assert "TIMEOUT_FAILURE" in fw_content
    # branch name derived from MiniPRD_HungBranch.md → feature/miniprd-hungbranch
    assert "miniprd-hungbranch" in fw_content


# ---------------------------------------------------------------------------
# Test 7 — SIGTERM handler releases lock and writes partial FAILED_WORKFLOWS.md
# ---------------------------------------------------------------------------

def test_sigterm_releases_lock(tmp_path, monkeypatch):
    """SIGTERM handler must release fcntl lock and write entries for in-flight branches."""
    monkeypatch.chdir(tmp_path)
    fw = _patch_fw(monkeypatch, tmp_path)
    monkeypatch.setattr(orch, "FAILED_WORKFLOWS_PATH", fw)
    monkeypatch.setattr(orch, "LOG_DIR", tmp_path / ".agents" / "logs")

    lock_path = tmp_path / ".orchestrator.lock"
    lock_fd = open(lock_path, "w")

    # Simulate in-flight state
    orch._lock_fd = lock_fd
    orch._in_flight_branches = ["feature/branch-a", "feature/branch-b"]
    orch._structured_log_path = tmp_path / ".agents" / "logs" / "orch.jsonl"
    orch._structured_log_path.parent.mkdir(parents=True, exist_ok=True)

    # Call cleanup handler directly (simulates SIGTERM)
    try:
        t0 = time.monotonic()
        with pytest.raises(SystemExit):
            orch._cleanup_handler(signal.SIGTERM, None)
        elapsed = time.monotonic() - t0
        assert elapsed < 5.0, f"SIGTERM handler took {elapsed:.2f}s (must complete within 5s)"
    finally:
        # Restore global state
        orch._in_flight_branches = []
        orch._lock_fd = None

    assert fw.exists()
    fw_content = fw.read_text()
    assert "feature/branch-a" in fw_content
    assert "feature/branch-b" in fw_content
    lock_fd.close()


# ---------------------------------------------------------------------------
# Test 8 — Idempotency: node with existing _provenance is skipped
# ---------------------------------------------------------------------------

def test_idempotency_skips_provisioned_node(tmp_path, monkeypatch):
    """
    If a node already has _provenance in architecture.yml, _worker must return
    SKIPPED_PROVENANCE_EXISTS without doing any work.
    """
    monkeypatch.chdir(tmp_path)
    arch_path = tmp_path / "architecture.yml"
    node = {
        "id": "miniprd_test",
        "dimension": "Module",
        "status": "clean",
        "_provenance": {
            "branch_name": "feature/existing",
            "timestamp": "2025-01-01T00:00:00Z",
        },
    }
    arch_path.write_text(
        yaml.dump({"nodes": [node]}, sort_keys=False), encoding="utf-8"
    )
    monkeypatch.setattr(orch, "ARCH_YAML_PATH", arch_path)

    miniprd = tmp_path / "MiniPRD_Test.md"
    miniprd.write_text("# MiniPRD\n", encoding="utf-8")

    # node_id derived from MiniPRD_Test.md → "miniprd_test"
    assert orch._has_provenance("miniprd_test", arch_path) is True
    assert orch._has_provenance("other_node", arch_path) is False


# ---------------------------------------------------------------------------
# Test 9 — Orchestrator-level TOKEN_OVERRUN when API returns stop_reason=max_tokens
# ---------------------------------------------------------------------------

def test_orchestrator_token_overrun_from_api(tmp_path, monkeypatch):
    """When Anthropic API returns stop_reason='max_tokens', _worker must write
    TOKEN_OVERRUN to FAILED_WORKFLOWS and NOT call the daemon subprocess."""
    monkeypatch.chdir(tmp_path)
    fw = _patch_fw(monkeypatch, tmp_path)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key")

    # Node WITHOUT _provenance so idempotency check passes through.
    # _node_id_for_miniprd("MiniPRD_Test.md") → strip "MiniPRD_" → lower → "test"
    arch_path = tmp_path / "architecture.yml"
    arch_path.write_text(
        yaml.dump(
            {"nodes": [{"id": "test", "dimension": "Module", "status": "clean"}]},
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    miniprd = tmp_path / "MiniPRD_Test.md"
    miniprd.write_text("# MiniPRD\n", encoding="utf-8")

    daemon_called = []

    def _fake_subprocess_run(cmd, **kwargs):
        if isinstance(cmd, list) and any("hyper_daemon" in str(c) for c in cmd):
            daemon_called.append(list(cmd))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(orch.subprocess, "run", _fake_subprocess_run)

    mock_message = SimpleNamespace(id="msg_test_overrun", stop_reason="max_tokens")
    mock_client = MagicMock()
    mock_client.messages.create.return_value = mock_message
    mock_anthropic_cls = MagicMock(return_value=mock_client)
    monkeypatch.setattr(orch.anthropic, "Anthropic", mock_anthropic_cls)

    monkeypatch.setattr(orch, "_resolve_model", lambda skill_name: "claude-sonnet-4-6")

    result = orch._worker(
        superprd_name="TestSuperPRD",
        superprd_content="# SuperPRD\n",
        miniprd_path=miniprd,
        miniprd_content="# MiniPRD\n",
        arch_yaml_path=arch_path,
    )

    assert daemon_called == [], f"Daemon must not be called on TOKEN_OVERRUN; got: {daemon_called}"
    assert fw.exists(), "FAILED_WORKFLOWS.md must exist"
    assert "TOKEN_OVERRUN" in fw.read_text()
    assert result["success"] is False
    assert result["error"] == "TOKEN_OVERRUN"


# ---------------------------------------------------------------------------
# Test 10 — Idempotency: _worker returns immediately when _provenance exists
# ---------------------------------------------------------------------------

def test_idempotency_worker_skips_e2e(tmp_path, monkeypatch):
    """When architecture.yml already has _provenance for the target node, _worker
    must return the skip result dict without invoking any subprocess."""
    monkeypatch.chdir(tmp_path)

    arch_path = tmp_path / "architecture.yml"
    node_with_provenance = {
        "id": "test",
        "dimension": "Module",
        "status": "clean",
        "_provenance": {
            "branch_name": "feature/miniprd-test",
            "timestamp": "2025-01-01T00:00:00+00:00",
        },
    }
    arch_path.write_text(
        yaml.dump({"nodes": [node_with_provenance]}, sort_keys=False),
        encoding="utf-8",
    )
    miniprd = tmp_path / "MiniPRD_Test.md"
    miniprd.write_text("# MiniPRD\n", encoding="utf-8")

    def _must_not_call(*args, **kwargs):
        raise AssertionError(
            f"subprocess.run must not be called when idempotency triggers; args={args}"
        )

    monkeypatch.setattr(orch.subprocess, "run", _must_not_call)

    result = orch._worker(
        superprd_name="test",
        superprd_content="# SuperPRD\n",
        miniprd_path=miniprd,
        miniprd_content="# MiniPRD\n",
        arch_yaml_path=arch_path,
    )

    assert result["skipped"] is True
    assert result["error"] == "SKIPPED_PROVENANCE_EXISTS"
    assert result["success"] is False
