#!/usr/bin/env python3
"""
hyper_update_core.py — Core logic for /hyper-update skill

Implements Phase A-E:
- Phase A: Upstream fetching & GPG verification
- Phase B: Mutex lock management
- Phase C: File change detection
- Phase D: Auto-update non-sensitive files
- Phase E: Preflight validation
"""

import os
import sys
import json
import tempfile
import shutil
import subprocess
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional

# Import verification module
from hyper_update_verification import UpstreamVerification


class HyperUpdateCore:
    """Core logic orchestrator for framework upgrades."""

    LOCK_FILE = ".agents/.hyper-update.lock"
    UPSTREAM_REPO = "https://github.com/tjmustard/Hypergraph-Coding-Agent-Framework.git"

    # Non-sensitive files that get auto-updated
    AUTO_UPDATE_DIRS = [
        ".agents/skills/",
        ".agents/scripts/",
        ".agents/schemas/",
        "spec/",
        "tests/",
    ]
    AUTO_UPDATE_FILES = [".agentignore"]

    # Sensitive files requiring user approval (replaced/merged/skipped)
    SENSITIVE_FILES = [
        "CLAUDE.md",
        "AGENTS.md",
        "GEMINI.md",
    ]

    # Maps local destination filename -> upstream source path.
    # These files are maintained as install templates so that the installed versions
    # are framed for user projects, while the repo root files govern HACF development.
    SENSITIVE_FILE_UPSTREAM_SOURCES = {
        "CLAUDE.md": ".agents/install-templates/CLAUDE.md",
        "AGENTS.md": ".agents/install-templates/AGENTS.md",
        "GEMINI.md": ".agents/install-templates/GEMINI.md",
    }
    SENSITIVE_DIRS = [
        ".claude/commands/",
        ".clinerules/",
        ".roo/rules/",
        ".roo/rules-code/",
        ".cursor/rules/",
        ".windsurf/rules/",
    ]

    def __init__(self, project_root: str = "."):
        self.project_root = Path(project_root).resolve()
        self.upstream_dir = None
        self.commit_sha = None
        self.gpg_verified = False
        self.backups_cleaned = 0
        self.log_file = None
        self.verifier = UpstreamVerification(str(self.project_root))
        self.changes = {
            "auto_updated": [],
            "replaced": [],
            "merged": [],
            "skipped": [],
            "files_with_changes": {},  # For later processing
        }

    def log(self, msg: str):
        """Log a message with timestamp."""
        print(f"[hyper-update] {msg}", file=sys.stderr)

    def error(self, msg: str):
        """Log an error and exit."""
        self.log(f"❌ Error: {msg}")
        sys.exit(1)

    def run_cmd(self, cmd: List[str], check: bool = True, capture: bool = False) -> str:
        """Run a shell command and return output."""
        try:
            result = subprocess.run(
                cmd,
                cwd=str(self.project_root),
                check=check,
                capture_output=capture,
                text=True,
            )
            return result.stdout.strip() if capture else ""
        except subprocess.CalledProcessError as e:
            if check:
                self.error(f"Command failed: {' '.join(cmd)}\n{e.stderr}")
            return ""

    # ==== PHASE A: Upstream Fetching & Verification ====

    def verify_git_available(self):
        """Verify git is installed."""
        self.log("Checking git availability...")
        try:
            subprocess.run(["git", "--version"], check=True, capture_output=True)
            self.log("✅ Git available")
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.error("git not found. Install git and retry.")

    def fetch_upstream(self):
        """Fetch upstream repo into temp directory."""
        self.log("Fetching upstream repository...")
        try:
            tmpdir = tempfile.mkdtemp(prefix="hcaf-upstream-")
            self.upstream_dir = Path(tmpdir)

            subprocess.run(
                [
                    "git", "clone",
                    "--depth=1",
                    "--branch", "main",
                    self.UPSTREAM_REPO,
                    str(self.upstream_dir),
                ],
                check=True,
                capture_output=True,
            )
            self.log(f"✅ Upstream fetched to {tmpdir}")
        except subprocess.CalledProcessError as e:
            self.error(f"Failed to fetch upstream: {e.stderr}")

    def verify_gpg_signature(self):
        """Verify GPG signature of upstream HEAD commit."""
        self.log("Verifying GPG signature...")

        # Check if GPG is available first
        if not self.verifier.check_gpg_available():
            self.error(
                "GPG is required but not installed.\n"
                "Install GPG with: apt-get install gnupg (Linux) / brew install gnupg (macOS)\n"
                "Abort."
            )

        # Get the commit SHA
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(self.upstream_dir),
            check=True,
            capture_output=True,
            text=True,
        )
        self.commit_sha = result.stdout.strip()
        self.log(f"Upstream commit: {self.commit_sha[:8]}")

        # Verify the signature
        result = subprocess.run(
            ["git", "verify-commit", self.commit_sha],
            cwd=str(self.upstream_dir),
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            self.gpg_verified = True
            self.log("✅ GPG signature verified")
        else:
            self.error("Upstream commit is not GPG-signed or signature is invalid.\n"
                      "Upstream repository requires signed commits for security.\n"
                      "Abort.")

    # ==== PHASE B: Mutex Lock Management ====

    def acquire_lock(self):
        """Create mutex lock file."""
        lock_path = self.project_root / self.LOCK_FILE

        if lock_path.exists():
            self.error(
                "Another upgrade in progress. Wait for it to complete or "
                "delete .agents/.hyper-update.lock manually."
            )

        lock_path.parent.mkdir(parents=True, exist_ok=True)
        lock_path.write_text(f"locked at {datetime.utcnow().isoformat()}\n")
        self.log("🔒 Lock acquired")

    def release_lock(self):
        """Delete mutex lock file."""
        lock_path = self.project_root / self.LOCK_FILE
        if lock_path.exists():
            lock_path.unlink()
            self.log("🔓 Lock released")

    # ==== PHASE C: File Change Detection ====

    def detect_changes(self) -> Dict[str, List[str]]:
        """Detect which files have changed in upstream."""
        changes = {
            "auto_update": [],
            "sensitive": [],
            "unchanged": [],
        }

        self.log("Detecting upstream changes...")

        # Check auto-update directories
        for auto_dir in self.AUTO_UPDATE_DIRS:
            upstream_path = self.upstream_dir / auto_dir
            local_path = self.project_root / auto_dir

            if upstream_path.exists():
                # Check if directory has any files (simple existence check)
                if list(upstream_path.rglob("*")):
                    changes["auto_update"].append(auto_dir)

        # Check auto-update files
        for auto_file in self.AUTO_UPDATE_FILES:
            upstream_path = self.upstream_dir / auto_file
            local_path = self.project_root / auto_file

            if upstream_path.exists() and local_path.exists():
                # Check if file has changed
                if not self._files_identical(upstream_path, local_path):
                    changes["auto_update"].append(auto_file)

        # Check sensitive files
        for sens_file in self.SENSITIVE_FILES:
            upstream_src = self.SENSITIVE_FILE_UPSTREAM_SOURCES.get(sens_file, sens_file)
            upstream_path = self.upstream_dir / upstream_src
            local_path = self.project_root / sens_file

            if upstream_path.exists() and local_path.exists():
                if not self._files_identical(upstream_path, local_path):
                    changes["sensitive"].append(sens_file)
                    self.changes["files_with_changes"][sens_file] = {
                        "upstream": str(upstream_path),
                        "local": str(local_path),
                    }

        # Check sensitive directories
        for sens_dir in self.SENSITIVE_DIRS:
            upstream_dir_path = self.upstream_dir / sens_dir
            local_dir_path = self.project_root / sens_dir

            if upstream_dir_path.exists() and local_dir_path.exists():
                for upstream_file in upstream_dir_path.rglob("*"):
                    if upstream_file.is_file():
                        rel_path = upstream_file.relative_to(self.upstream_dir)
                        local_file = self.project_root / rel_path

                        if local_file.exists():
                            if not self._files_identical(upstream_file, local_file):
                                changes["sensitive"].append(str(rel_path))
                                self.changes["files_with_changes"][str(rel_path)] = {
                                    "upstream": str(upstream_file),
                                    "local": str(local_file),
                                }

        return changes

    def _files_identical(self, file1: Path, file2: Path) -> bool:
        """Check if two files are identical."""
        try:
            return file1.read_bytes() == file2.read_bytes()
        except (OSError, IOError):
            return False

    # ==== PHASE D: Auto-Update Non-Sensitive Files ====

    def auto_update_files(self):
        """Auto-update non-sensitive files without prompting."""
        self.log("Auto-updating non-sensitive files...")

        for auto_dir in self.AUTO_UPDATE_DIRS:
            upstream_path = self.upstream_dir / auto_dir
            local_path = self.project_root / auto_dir

            if upstream_path.exists():
                if local_path.exists():
                    shutil.rmtree(local_path)
                shutil.copytree(upstream_path, local_path)
                self.changes["auto_updated"].append(auto_dir)
                self.log(f"✅ Auto-updated: {auto_dir}")

        for auto_file in self.AUTO_UPDATE_FILES:
            upstream_path = self.upstream_dir / auto_file
            local_path = self.project_root / auto_file

            if upstream_path.exists():
                shutil.copy2(upstream_path, local_path)
                self.changes["auto_updated"].append(auto_file)
                self.log(f"✅ Auto-updated: {auto_file}")

    # ==== PHASE E: Preflight Validation ====

    def check_uncommitted_changes(self) -> List[str]:
        """Check for uncommitted changes in sensitive areas. Returns list of affected files."""
        self.log("Checking for uncommitted changes...")

        result = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=str(self.project_root),
            check=False,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            self.log("⚠️  Could not run git status (repo may not be initialized)")
            return []

        # Parse uncommitted changes
        uncommitted = []
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            status, filepath = line[:2], line[3:]

            # Check if it's in sensitive areas
            if any(
                filepath.startswith(path)
                for path in ["CLAUDE.md", "AGENTS.md", "GEMINI.md", ".agents/", ".claude/"]
            ):
                uncommitted.append(filepath)

        return uncommitted

    # ==== Main Execution ====

    def run(self) -> Dict:
        """Execute the full upgrade sequence."""
        try:
            # Phase A: Upstream Fetching & Verification
            self.verify_git_available()
            self.fetch_upstream()
            self.verify_gpg_signature()

            # Phase: Cleanup old backups (part of upstream verification)
            self.backups_cleaned, cleanup_msg = self.verifier.cleanup_old_backups()
            self.log(cleanup_msg)

            # Phase: Get backup age warnings
            warnings = self.verifier.get_backup_age_warnings()

            # Phase B: Acquire lock
            self.acquire_lock()

            # Phase E: Preflight validation (check for uncommitted changes)
            uncommitted = self.check_uncommitted_changes()
            if uncommitted:
                # Log the warning and return early; let skill layer handle user interaction
                self.log_file = self.verifier.log_audit_trail(
                    status="preflight_warning",
                    commit_sha=self.commit_sha,
                    gpg_signature="verified",
                    backups_cleaned=self.backups_cleaned,
                    error_msg="uncommitted changes detected",
                )
                result = {
                    "status": "preflight_warning",
                    "commit_sha": self.commit_sha,
                    "gpg_verified": self.gpg_verified,
                    "uncommitted_changes": uncommitted,
                    "upstream_dir": str(self.upstream_dir),
                    "backups_cleaned": self.backups_cleaned,
                    "log_file": self.log_file,
                    "backup_warnings": warnings,
                }
                print(json.dumps(result, indent=2))
                return result

            # Phase C: Detect changes
            changes = self.detect_changes()
            self.log(f"Found {len(changes['sensitive'])} sensitive files with changes")

            # Phase D: Auto-update non-sensitive files
            if changes["auto_update"]:
                self.auto_update_files()
            else:
                self.log("No non-sensitive files to update")

            # Log successful completion
            total_files = len(self.changes["auto_updated"]) + len(self.changes["replaced"]) + len(self.changes["merged"])
            self.log_file = self.verifier.log_audit_trail(
                status="success",
                commit_sha=self.commit_sha,
                gpg_signature="verified",
                backups_cleaned=self.backups_cleaned,
                files_processed=total_files,
            )

            # Return results for next phases
            result = {
                "status": "success",
                "commit_sha": self.commit_sha,
                "gpg_verified": self.gpg_verified,
                "changes": self.changes,
                "sensitive_files": changes["sensitive"],
                "upstream_dir": str(self.upstream_dir),
                "backups_cleaned": self.backups_cleaned,
                "log_file": self.log_file,
                "backup_warnings": warnings,
            }

            print(json.dumps(result, indent=2))
            return result

        except Exception as e:
            self.error(f"Unexpected error: {str(e)}")
        finally:
            # Always clean up lock and temp directory
            self.release_lock()
            if self.upstream_dir and self.upstream_dir.exists():
                shutil.rmtree(self.upstream_dir)


def main():
    """Entry point."""
    project_root = os.environ.get("PROJECT_ROOT", ".")
    updater = HyperUpdateCore(project_root)
    updater.run()


if __name__ == "__main__":
    main()
