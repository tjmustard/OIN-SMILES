#!/usr/bin/env python3
"""
Hyper Daemon — Execute → Audit → Oracle → Fix loop

Supports two modes:
  - Standalone:     hyper_daemon.py <superprd_name> <miniprd_path>
  - Orchestrated:   hyper_daemon.py --orchestrated <superprd_name> <miniprd_path> [--api-response-id <id>]

The standalone mode behavior is identical in both paths (same loop, same guards).
"""

import sys
if sys.version_info < (3, 11):
    sys.exit("Python 3.11+ required")

import argparse
import ast
import difflib
import os
import re
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MAX_ITERATIONS = 5
FAILED_WORKFLOWS_PATH = Path("spec/active/FAILED_WORKFLOWS.md")
FAILED_WORKFLOWS_ARCHIVE_DIR = Path("spec/archive")
SPEC_DRIFT_THRESHOLD = 0.40  # 40% drift triggers SPEC_DRIFT

# Sentinel exit code emitted by hyper_fix.py on TOKEN_OVERRUN; daemon checks
# this to break the loop without counting the aborted call as an iteration
EXIT_TOKEN_OVERRUN = 3


# ---------------------------------------------------------------------------
# Failure type enum (mirrors hyper_orchestrator.py — shared truth via enum)
# ---------------------------------------------------------------------------

class FailureType(str, Enum):
    ORACLE_FAILURE = "ORACLE_FAILURE"
    TOKEN_OVERRUN = "TOKEN_OVERRUN"
    NO_TESTS_COLLECTED = "NO_TESTS_COLLECTED"
    SPEC_DRIFT = "SPEC_DRIFT"
    UNSUPPORTED_CONFLICT_TYPE = "UNSUPPORTED_CONFLICT_TYPE"
    RESOLVER_MALFORMED_RESPONSE = "RESOLVER_MALFORMED_RESPONSE"
    CONTEXT_OVERFLOW = "CONTEXT_OVERFLOW"
    TIMEOUT_FAILURE = "TIMEOUT_FAILURE"
    REBASE_CONFLICT = "REBASE_CONFLICT"


# ---------------------------------------------------------------------------
# FAILED_WORKFLOWS.md helpers
# ---------------------------------------------------------------------------

def _count_failed_headers(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        return text.count("## [FAILED]")
    except OSError:
        return 0


def _rotate_failed_workflows() -> None:
    """Archive FAILED_WORKFLOWS.md when it reaches 50 entries."""
    import shutil
    FAILED_WORKFLOWS_ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    iso = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = FAILED_WORKFLOWS_ARCHIVE_DIR / f"FAILED_WORKFLOWS_{iso}.md"
    shutil.copy2(FAILED_WORKFLOWS_PATH, dest)
    FAILED_WORKFLOWS_PATH.write_text("# Failed Workflows\n\n", encoding="utf-8")
    print(f"[ROTATE] FAILED_WORKFLOWS.md archived to {dest}")


def append_failed_workflow(
    branch_name: str,
    superprd_name: str,
    miniprd_path: str,
    iterations_attempted: int,
    failure_type: FailureType,
    oracle_exit_code: int,
    detail: str = "",
) -> None:
    """Append a structured failure entry; rotate at 50 entries."""
    FAILED_WORKFLOWS_PATH.parent.mkdir(parents=True, exist_ok=True)
    if not FAILED_WORKFLOWS_PATH.exists():
        FAILED_WORKFLOWS_PATH.write_text("# Failed Workflows\n\n", encoding="utf-8")

    if _count_failed_headers(FAILED_WORKFLOWS_PATH) >= 50:
        _rotate_failed_workflows()

    iso = datetime.now(timezone.utc).isoformat()
    entry = (
        f"\n## [FAILED] {branch_name} — {iso}\n\n"
        f"- **SuperPRD:** {superprd_name}\n"
        f"- **MiniPRD:** {miniprd_path}\n"
        f"- **Branch:** {branch_name}\n"
        f"- **Iterations Attempted:** {iterations_attempted}/5\n"
        f"- **Failure Type:** {failure_type.value}\n"
        f"- **Final Oracle Exit Code:** {oracle_exit_code}\n"
    )
    if detail:
        entry += f"\n{detail}\n"
    entry += "\n---\n"

    with FAILED_WORKFLOWS_PATH.open("a", encoding="utf-8") as fh:
        fh.write(entry)


# ---------------------------------------------------------------------------
# Agent subprocess runner
# ---------------------------------------------------------------------------

def _run_agent_script(
    script_name: str,
    args: list[str],
    env_vars: dict | None = None,
) -> bool:
    """Execute an existing framework script. Returns True on returncode==0."""
    env = os.environ.copy()
    if env_vars:
        env.update(env_vars)
    script_path = Path(__file__).parent / f"{script_name}.py"
    cmd = [sys.executable, str(script_path)] + args
    print(f"    [AGENT] Executing: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, env=env, capture_output=False, check=False)
    if result.returncode != 0:
        print(f"    [AGENT] {script_name} exited with code {result.returncode}")
    return result.returncode == 0


# ---------------------------------------------------------------------------
# Task 18 — Oracle: --collect-only pre-flight
# ---------------------------------------------------------------------------

def _pytest_collect_count() -> int:
    """
    Run pytest --collect-only -q and parse item count.
    Returns the number of collected test items.
    """
    result = subprocess.run(
        ["uv", "run", "pytest", "--collect-only", "-q"],
        capture_output=True, text=True, check=False,
    )
    output = result.stdout + result.stderr
    # pytest outputs lines like "3 tests collected" or "no tests ran"
    match = re.search(r"(\d+)\s+(?:test[s]?\s+)?(?:collected|item[s]?)", output, re.IGNORECASE)
    if match:
        return int(match.group(1))
    # If output says "no tests ran" or "0 tests" explicitly
    if re.search(r"\bno tests\b|\b0 test|\b0 item", output, re.IGNORECASE):
        return 0
    # If no match and exit code is 5 (no tests collected by pytest convention)
    if result.returncode == 5:
        return 0
    # Default: assume non-zero if we can't parse (safer than blocking)
    return 1


def _run_oracle_tests() -> tuple[int, str, str]:
    """
    Run pytest with --maxfail=1 --tb=short.
    Returns (exit_code, stdout, stderr).
    """
    result = subprocess.run(
        ["uv", "run", "pytest", "--maxfail=1", "--tb=short"],
        capture_output=True, text=True, check=False,
    )
    return result.returncode, result.stdout, result.stderr


# ---------------------------------------------------------------------------
# Task 19 — Spec drift detector
# ---------------------------------------------------------------------------

def _extract_public_interface(source: str) -> set[str]:
    """
    Extract top-level public function/class names from Python source.
    Checks __all__ first; falls back to top-level def/class with no leading underscore.
    """
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return set()

    # Check for explicit __all__
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__all__":
                    if isinstance(node.value, (ast.List, ast.Tuple)):
                        names: set[str] = set()
                        for elt in node.value.elts:
                            if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                                names.add(elt.value)
                        return names

    # Fall back: top-level defs and classes that don't start with _
    names = set()
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if not node.name.startswith("_"):
                names.add(node.name)
    return names


def _check_spec_drift(
    miniprd_path: Path,
    baseline_interface: set[str],
) -> bool:
    """
    Return True if spec drift exceeds SPEC_DRIFT_THRESHOLD.
    Compares current on-disk Python files against the baseline interface captured
    at iteration 1.
    """
    if not baseline_interface:
        return False

    # Find Python files that the MiniPRD might generate (same directory as outputs)
    # Heuristic: scan all .py files touched since daemon start in working directory
    current_names: set[str] = set()
    try:
        result = subprocess.run(
            ["git", "diff", "--name-only", "HEAD"],
            capture_output=True, text=True, check=False,
        )
        changed_files = [
            f for f in result.stdout.splitlines()
            if f.endswith(".py")
        ]
        for filepath in changed_files:
            p = Path(filepath)
            if p.exists():
                try:
                    src = p.read_text(encoding="utf-8", errors="replace")
                    current_names |= _extract_public_interface(src)
                except OSError:
                    pass
    except Exception:
        return False

    if not current_names:
        return False

    # MiniPRD specifies "diff-line similarity" metric (difflib.SequenceMatcher)
    baseline_text = "\n".join(sorted(baseline_interface))
    current_text = "\n".join(sorted(current_names))
    ratio = difflib.SequenceMatcher(None, baseline_text, current_text).ratio()
    return (1 - ratio) > SPEC_DRIFT_THRESHOLD


def _capture_baseline_interface(miniprd_path: Path) -> set[str]:
    """
    Capture the public interface from any Python files modified after the first execute.
    Called at the end of iteration 1.
    """
    names: set[str] = set()
    try:
        result = subprocess.run(
            ["git", "diff", "--name-only", "HEAD"],
            capture_output=True, text=True, check=False,
        )
        for filepath in result.stdout.splitlines():
            if filepath.endswith(".py"):
                p = Path(filepath)
                if p.exists():
                    try:
                        src = p.read_text(encoding="utf-8", errors="replace")
                        names |= _extract_public_interface(src)
                    except OSError:
                        pass
    except Exception:
        pass
    return names


# ---------------------------------------------------------------------------
# Core Execute → Audit → Oracle → Fix loop
# ---------------------------------------------------------------------------

def run_loop(
    superprd_name: str,
    miniprd_path: Path,
    orchestrated_mode: bool = False,
    api_response_id: str = "",
) -> bool:
    """
    Execute the full autonomous loop.
    Returns True on success (all tests pass within MAX_ITERATIONS).

    Tasks implemented here:
      - Task 17: stop_reason guard (orchestrated mode only — API already called upstream)
      - Task 18: --collect-only pre-flight before each pytest call
      - Task 19: spec drift detection on iterations 2–5
      - Task 20: iteration hard cap at 5

    In orchestrated mode the Anthropic API call was already made by the orchestrator;
    the daemon loop focuses on Execute → Audit → Oracle → Fix driven by subprocess calls.
    """
    branch_name = _derive_branch_name(miniprd_path)
    env_vars = {"HYPERGRAPH_PROVENANCE": superprd_name}
    if orchestrated_mode and api_response_id:
        env_vars["HYPER_API_RESPONSE_ID"] = api_response_id

    baseline_interface: set[str] = set()
    iteration = 0
    oracle_exit_code = 1  # tracks last oracle result for ORACLE_FAILURE delegation

    # Task 20 — iteration hard cap: loop body only runs while iteration < MAX_ITERATIONS
    while iteration < MAX_ITERATIONS:
        iteration += 1
        print(f"\n[DAEMON] Iteration {iteration}/{MAX_ITERATIONS}: {miniprd_path.name}")

        # --- Task 19: Spec drift check at START of each iteration >= 2 (before execute) ---
        if iteration >= 2 and baseline_interface:
            if _check_spec_drift(miniprd_path, baseline_interface):
                print("[DAEMON] SPEC_DRIFT detected: >40% diff-line similarity change.")
                append_failed_workflow(
                    branch_name=branch_name,
                    superprd_name=superprd_name,
                    miniprd_path=str(miniprd_path),
                    iterations_attempted=iteration - 1,  # this iteration not consumed
                    failure_type=FailureType.SPEC_DRIFT,
                    oracle_exit_code=oracle_exit_code,
                    detail=(
                        f"More than {int(SPEC_DRIFT_THRESHOLD * 100)}% diff-line similarity "
                        "change detected relative to the MiniPRD spec baseline. "
                        "Manual review required."
                    ),
                )
                return False

        # --- Step 1: Execute (first iteration only; fix agent patches in-place after) ---
        if iteration == 1:
            execute_ok = _run_agent_script(
                "hyper_execute_runner",
                [str(miniprd_path)],
                env_vars,
            )
            if not execute_ok:
                print(f"[DAEMON] Execute agent failed on iteration {iteration}. Retrying…")
                continue

            # Capture baseline interface after first execute (for drift detection)
            baseline_interface = _capture_baseline_interface(miniprd_path)

        # --- Step 2: Audit ---
        audit_ok = _run_agent_script(
            "hyper_audit_runner",
            [str(miniprd_path)],
            env_vars,
        )
        if not audit_ok:
            print(f"[DAEMON] Audit failed on iteration {iteration}. Engaging fix…")
            # Fall through to oracle anyway; fix will handle both audit and test failures

        # --- Task 18: --collect-only pre-flight ---
        test_count = _pytest_collect_count()
        if test_count == 0:
            print("[DAEMON] pytest --collect-only found 0 tests. Aborting.")
            append_failed_workflow(
                branch_name=branch_name,
                superprd_name=superprd_name,
                miniprd_path=str(miniprd_path),
                iterations_attempted=iteration - 1,  # do NOT count this as an iteration
                failure_type=FailureType.NO_TESTS_COLLECTED,
                oracle_exit_code=5,
                detail="pytest --collect-only reported 0 test items. "
                       "Verify test files exist and are importable.",
            )
            return False

        # --- Step 3: Oracle ---
        oracle_exit_code, oracle_stdout, oracle_stderr = _run_oracle_tests()

        if oracle_exit_code == 0:
            print(f"[DAEMON] All tests passed on iteration {iteration}!")
            return True

        print(f"[DAEMON] Oracle failed (exit {oracle_exit_code}) on iteration {iteration}.")

        # --- Step 4: Fix (only if more iterations remain) ---
        if iteration < MAX_ITERATIONS:
            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".log",
                prefix=f"oracle_err_{iteration}_",
                delete=False,
                encoding="utf-8",
            ) as tmp_log:
                tmp_log.write(
                    f"EXIT CODE: {oracle_exit_code}\n"
                    f"STDOUT:\n{oracle_stdout}\n"
                    f"STDERR:\n{oracle_stderr}\n"
                )
                tmp_log_path = tmp_log.name

            try:
                env = os.environ.copy()
                env.update(env_vars)
                fix_script = Path(__file__).parent / "hyper_fix.py"
                fix_result = subprocess.run(
                    [
                        sys.executable, str(fix_script),
                        str(miniprd_path),
                        "--error-log", tmp_log_path,
                        "--branch", branch_name,
                        "--iteration", str(iteration),
                    ],
                    env=env, capture_output=False, check=False,
                )
                # Task 17 — stop_reason guard (daemon implementation):
                # The daemon never calls client.messages.create() directly; API calls happen
                # in hyper_orchestrator._worker() (orchestrated mode) and inside hyper_fix.py
                # (both modes). Both upstream callers detect stop_reason=="max_tokens" and
                # exit with EXIT_TOKEN_OVERRUN (3). This sentinel catch is the daemon's Task
                # 17 contract: it enforces no-file-write and FAILED_WORKFLOWS logging for any
                # TOKEN_OVERRUN that propagates from any upstream source.
                if fix_result.returncode == EXIT_TOKEN_OVERRUN:
                    print("[DAEMON] TOKEN_OVERRUN from fix agent — aborting loop.")
                    append_failed_workflow(
                        branch_name=branch_name,
                        superprd_name=superprd_name,
                        miniprd_path=str(miniprd_path),
                        iterations_attempted=iteration,
                        failure_type=FailureType.TOKEN_OVERRUN,
                        oracle_exit_code=oracle_exit_code,
                        detail="hyper_fix exited EXIT_TOKEN_OVERRUN (stop_reason=max_tokens); "
                               "truncated output not written to disk.",
                    )
                    return False
            finally:
                try:
                    os.unlink(tmp_log_path)
                except OSError:
                    pass

    # Task 20 — max iterations reached; delegate ORACLE_FAILURE write to hyper_fix.py
    # so it is the canonical owner (Task 6 Autonomous Resolution MiniPRD)
    print(f"[DAEMON] ORACLE_FAILURE: max {MAX_ITERATIONS} iterations reached without passing tests.")
    hyper_fix_script = Path(__file__).parent / "hyper_fix.py"
    env = os.environ.copy()
    env.update(env_vars)
    subprocess.run(
        [
            sys.executable, str(hyper_fix_script),
            str(miniprd_path), "--error-log", os.devnull,
            "--branch", branch_name,
            "--oracle-exit-code", str(oracle_exit_code),
            "--record-oracle-failure",
        ],
        env=env, check=False,
    )
    return False


# ---------------------------------------------------------------------------
# Orchestrated mode — called by hyper_orchestrator (Task 16)
# ---------------------------------------------------------------------------

def run_orchestrated_loop(
    superprd_name: str,
    miniprd_path: Path,
    api_response_id: str = "",
) -> bool:
    """
    Run the Execute → Audit → Oracle → Fix loop in orchestrated mode.
    The Anthropic API call was already made by the orchestrator; this function
    drives the subprocess-based execution pipeline only.
    """
    print(f"[DAEMON/ORCHESTRATED] Starting loop for {miniprd_path.name}")
    return run_loop(
        superprd_name=superprd_name,
        miniprd_path=miniprd_path,
        orchestrated_mode=True,
        api_response_id=api_response_id,
    )


# ---------------------------------------------------------------------------
# Standalone mode (Task 16)
# ---------------------------------------------------------------------------

def run_standalone(superprd_name: str, miniprd_path: Path) -> bool:
    """
    Standalone entry point. Behavior identical to orchestrated loop but
    without expecting a pre-existing API response ID.
    """
    print(f"[DAEMON/STANDALONE] Starting loop for {miniprd_path.name}")
    return run_loop(
        superprd_name=superprd_name,
        miniprd_path=miniprd_path,
        orchestrated_mode=False,
        api_response_id="",
    )


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _derive_branch_name(miniprd_path: Path) -> str:
    stem = miniprd_path.stem
    if stem.startswith("MiniPRD_"):
        stem = stem[len("MiniPRD_"):]
    slug = stem.lower().replace("_", "-").replace(" ", "-")
    return f"feature/{slug}"


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Hyper Daemon: Execute → Audit → Oracle → Fix loop."
    )
    parser.add_argument(
        "--orchestrated",
        action="store_true",
        default=False,
        help="Run in orchestrated mode (called by hyper_orchestrator).",
    )
    parser.add_argument(
        "superprd_name",
        type=str,
        help="SuperPRD name (used for provenance and FAILED_WORKFLOWS entries).",
    )
    parser.add_argument(
        "miniprd_path",
        type=Path,
        help="Path to the MiniPRD .md file to implement.",
    )
    parser.add_argument(
        "--api-response-id",
        type=str,
        default="",
        dest="api_response_id",
        help="Anthropic message ID from the orchestrator's API call (orchestrated mode only).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()

    # Task 16 — conditional branch entry point
    if args.orchestrated:
        success = run_orchestrated_loop(
            superprd_name=args.superprd_name,
            miniprd_path=args.miniprd_path,
            api_response_id=args.api_response_id,
        )
    else:
        success = run_standalone(
            superprd_name=args.superprd_name,
            miniprd_path=args.miniprd_path,
        )

    sys.exit(0 if success else 1)
