"""
Integration tests for MiniPRD: Autonomous Resolution

Covers:
  - hyper_resolve_conflict.py: whitelist, binary, JSON decode, confidence assertion,
    confidence gate, syntax override, full pass, error field
  - hyper_fix.py: spec drift
"""

import ast
import json
import os
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent.parent / ".agents" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import hyper_resolve_conflict as hrc  # noqa: E402
import hyper_fix as hf  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fake_message(content_text, stop_reason="end_turn"):
    msg = SimpleNamespace()
    msg.stop_reason = stop_reason
    msg.content = [SimpleNamespace(text=content_text)]
    return msg


def _resolve_response(confidence, resolved_code, error=None):
    payload = {"confidence_score": confidence, "resolved_code": resolved_code, "error": error}
    return json.dumps(payload)


# ---------------------------------------------------------------------------
# Test 1 — Whitelist blocks .lock file (no API call)
# ---------------------------------------------------------------------------

def test_whitelist_blocks_lock_file(tmp_path, monkeypatch):
    """A .lock conflict file must exit 1 with UNSUPPORTED_CONFLICT_TYPE, no API call."""
    lock_file = tmp_path / "poetry.lock"
    lock_file.write_text("[[package]]\nname = 'foo'\n", encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    api_called = []
    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.side_effect = (
            lambda *a, **kw: api_called.append(True) or _fake_message("{}")
        )
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(lock_file, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert not api_called, "API must NOT be called for unsupported file types"
    assert fw_path.exists()
    assert "UNSUPPORTED_CONFLICT_TYPE" in fw_path.read_text()


# ---------------------------------------------------------------------------
# Test 2 — Binary file blocked regardless of extension
# ---------------------------------------------------------------------------

def test_binary_file_blocked(tmp_path, monkeypatch):
    """A .py file containing null bytes must be treated as UNSUPPORTED_CONFLICT_TYPE."""
    binary_py = tmp_path / "module.py"
    # \x80\x81\x82 are not valid UTF-8 bytes (invalid leading bytes)
    binary_py.write_bytes(b"def foo():\n    pass\x80\x81\x82")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    with pytest.raises(SystemExit) as exc_info:
        hrc.resolve_conflict(binary_py, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert fw_path.exists()
    assert "UNSUPPORTED_CONFLICT_TYPE" in fw_path.read_text()


# ---------------------------------------------------------------------------
# Test 3 — JSON decode failure → RESOLVER_MALFORMED_RESPONSE
# ---------------------------------------------------------------------------

def test_json_decode_failure(tmp_path, monkeypatch):
    """Malformed JSON from API must set confidence=0, leave file unmodified, log RESOLVER_MALFORMED_RESPONSE."""
    py_file = tmp_path / "conflict.py"
    original = "def foo(): pass\n"
    py_file.write_text(original, encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.return_value = _fake_message("not json")
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(py_file, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert py_file.read_text() == original, "File must be unmodified"
    assert fw_path.exists()
    assert "RESOLVER_MALFORMED_RESPONSE" in fw_path.read_text()


# ---------------------------------------------------------------------------
# Test 4 — Confidence assertion: float value rejected
# ---------------------------------------------------------------------------

def test_confidence_assertion_float(tmp_path, monkeypatch):
    """confidence_score: 0.91 (float) must fail isinstance check → exit 1."""
    py_file = tmp_path / "f.py"
    py_file.write_text("x = 1\n", encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    response_json = json.dumps({"confidence_score": 0.91, "resolved_code": "x = 1\n", "error": None})
    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.return_value = _fake_message(response_json)
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(py_file, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert py_file.read_text() == "x = 1\n"
    assert fw_path.exists()
    assert "RESOLVER_MALFORMED_RESPONSE" in fw_path.read_text()


# ---------------------------------------------------------------------------
# Test 5 — Confidence assertion: out-of-range value rejected
# ---------------------------------------------------------------------------

def test_confidence_assertion_out_of_range(tmp_path, monkeypatch):
    """confidence_score: 150 must fail 0 <= score <= 100 check → exit 1."""
    py_file = tmp_path / "g.py"
    py_file.write_text("y = 2\n", encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    response_json = json.dumps({"confidence_score": 150, "resolved_code": "y = 2\n", "error": None})
    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.return_value = _fake_message(response_json)
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(py_file, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert py_file.read_text() == "y = 2\n"


# ---------------------------------------------------------------------------
# Test 6 — Confidence below gate (72 < 85) → file unmodified
# ---------------------------------------------------------------------------

def test_confidence_below_gate(tmp_path, monkeypatch):
    py_file = tmp_path / "h.py"
    original = "z = 3\n"
    py_file.write_text(original, encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    response_json = json.dumps({
        "confidence_score": 72,
        "resolved_code": "z = 999\n",
        "error": None,
        "reasoning": "low confidence explanation",
    })
    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.return_value = _fake_message(response_json)
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(py_file, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert py_file.read_text() == original


# ---------------------------------------------------------------------------
# Test 7 — Syntax check overrides high confidence
# ---------------------------------------------------------------------------

def test_syntax_overrides_confidence(tmp_path, monkeypatch):
    """confidence=92, but resolved_code has a SyntaxError → file unmodified."""
    py_file = tmp_path / "i.py"
    original = "a = 1\n"
    py_file.write_text(original, encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    bad_python = "def foo(:\n    pass\n"  # SyntaxError
    response_json = json.dumps({"confidence_score": 92, "resolved_code": bad_python, "error": None})
    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.return_value = _fake_message(response_json)
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(py_file, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert py_file.read_text() == original, "File must remain unmodified when syntax check fails"


# ---------------------------------------------------------------------------
# Test 8 — Full pass: confidence 90 + valid Python → file written, exit 0
# ---------------------------------------------------------------------------

def test_full_pass(tmp_path, monkeypatch):
    py_file = tmp_path / "j.py"
    py_file.write_text("<<<<<<< HEAD\nx = 1\n=======\nx = 2\n>>>>>>> branch\n", encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    resolved = "x = 2\n"
    response_json = json.dumps({"confidence_score": 90, "resolved_code": resolved, "error": None})
    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.return_value = _fake_message(response_json)
        # resolve_conflict calls sys.exit(0) on success — that's the expected exit code
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(py_file, branch="test", superprd_name="test", miniprd_path="test")
    assert exc_info.value.code == 0
    assert py_file.read_text() == resolved


# ---------------------------------------------------------------------------
# Test 9 — Spec drift in hyper_fix.py exits without consuming iteration
# ---------------------------------------------------------------------------

def test_spec_drift_exits_without_iteration(tmp_path, monkeypatch):
    """
    When >40% of MiniPRD-specified symbols are missing from current files,
    hyper_fix raises TerminalFailure with FT_SPEC_DRIFT, passing iteration-1
    so the counter is not incremented.
    """
    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hf, "FAILED_WORKFLOWS_PATH", fw_path)

    # MiniPRD with 5 code block symbols
    miniprd = tmp_path / "MiniPRD_Test.md"
    miniprd.write_text(
        "```python\ndef alpha(): pass\ndef beta(): pass\ndef gamma(): pass\n"
        "def delta(): pass\ndef epsilon(): pass\n```\n",
        encoding="utf-8",
    )

    # Current file only has 1 of the 5 symbols (80% missing > 40% threshold)
    cur_file = tmp_path / "cur.py"
    cur_file.write_text("def alpha(): pass\n", encoding="utf-8")

    error_log = tmp_path / "err.log"
    error_log.write_text("some error\n", encoding="utf-8")

    # Patch _collect_current_symbols to return only {alpha}
    monkeypatch.setattr(
        hf,
        "_collect_current_symbols",
        lambda content, path: {"alpha"},
    )

    with pytest.raises(hf.TerminalFailure) as exc_info:
        hf.run_fix(
            miniprd_path=str(miniprd),
            error_log_path=str(error_log),
            superprd_path=None,
            branch="feature/test",
            iteration=2,  # drift only checked on iterations >= 2
        )

    assert exc_info.value.failure_type == hf.FT_SPEC_DRIFT
    assert fw_path.exists()
    content = fw_path.read_text()
    assert "SPEC_DRIFT" in content
    # iterations_attempted should be iteration-1 = 1 (not consuming iteration 2)
    assert "1/5" in content


# ---------------------------------------------------------------------------
# Test 10 — Error field in response → RESOLVER_MALFORMED_RESPONSE
# ---------------------------------------------------------------------------

def test_error_field_resolver_malformed(tmp_path, monkeypatch):
    """response["error"] = "Content policy violation" → RESOLVER_MALFORMED_RESPONSE, exit 1."""
    py_file = tmp_path / "k.py"
    original = "pass\n"
    py_file.write_text(original, encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    fw_path = tmp_path / "spec" / "active" / "FAILED_WORKFLOWS.md"
    monkeypatch.setattr(hrc, "FAILED_WORKFLOWS_PATH", fw_path)

    response_json = json.dumps({
        "confidence_score": 90,
        "resolved_code": "pass\n",
        "error": "Content policy violation",
    })
    with patch("anthropic.Anthropic") as mock_anthropic:
        mock_anthropic.return_value.messages.create.return_value = _fake_message(response_json)
        with pytest.raises(SystemExit) as exc_info:
            hrc.resolve_conflict(py_file, branch="test", superprd_name="test", miniprd_path="test")

    assert exc_info.value.code == 1
    assert py_file.read_text() == original
    assert fw_path.exists()
    assert "RESOLVER_MALFORMED_RESPONSE" in fw_path.read_text()
