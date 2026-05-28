#!/usr/bin/env python3
"""
hyper_update_verification.py — Upstream verification, backup cleanup, and audit logging

Implements:
- GPG binary availability check
- Backup cleanup policy (30-day retention)
- Audit trail logging with JSON format
- Log rotation (keep last 10 log files)
"""

import os
import sys
import json
import shutil
import subprocess
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, Optional, Tuple


class UpstreamVerification:
    """Handles GPG verification, backup cleanup, and audit logging."""

    BACKUP_DIR = ".agents/.backup"
    LOG_DIR = ".agents/logs"
    LOG_FILE = "hyper-update.log"
    LOG_MAX_SIZE = 1024 * 1024  # 1MB
    LOG_MAX_FILES = 10
    BACKUP_RETENTION_DAYS = 30

    def __init__(self, project_root: str = "."):
        self.project_root = Path(project_root).resolve()
        self.backup_dir = self.project_root / self.BACKUP_DIR
        self.log_dir = self.project_root / self.LOG_DIR

    def log_message(self, msg: str):
        """Print a message to stderr for real-time feedback."""
        print(f"[hyper-update] {msg}", file=sys.stderr)

    def error(self, msg: str):
        """Log an error and exit."""
        self.log_message(f"❌ Error: {msg}")
        sys.exit(1)

    # ==== Phase: GPG Binary Check ====

    def check_gpg_available(self) -> bool:
        """
        Verify GPG is installed before attempting signature verification.
        Returns True if GPG is available, False otherwise.
        """
        self.log_message("Checking GPG availability...")
        try:
            result = subprocess.run(
                ["gpg", "--version"],
                check=True,
                capture_output=True,
                text=True,
            )
            self.log_message("✅ GPG is available")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def verify_upstream_commit(self, commit_sha: str, upstream_dir: Optional[str] = None) -> Tuple[bool, str]:
        """
        Verify GPG signature on upstream commit.
        Args: commit_sha (str), upstream_dir (optional path to cloned repo, defaults to current dir)
        Returns: (is_valid: bool, message: str)
        """
        cwd = upstream_dir if upstream_dir else str(self.project_root)

        try:
            result = subprocess.run(
                ["git", "verify-commit", commit_sha],
                cwd=cwd,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0:
                # Signature is valid
                sha_short = commit_sha[:7]
                return True, f"✅ GPG signature verified (commit SHA: {sha_short}...)"
            else:
                # Signature is invalid or missing
                return False, f"❌ Error: Upstream commit {commit_sha[:7]}... is not GPG-signed or signature is invalid.\nUpstream repository requires signed commits for security.\nAbort."
        except FileNotFoundError:
            return False, "❌ Error: git command not found"
        except Exception as e:
            return False, f"❌ Error: GPG verification failed: {e}"

    # ==== Phase: Backup Cleanup Policy ====

    def cleanup_old_backups(self) -> Tuple[int, str]:
        """
        Remove backup directories older than BACKUP_RETENTION_DAYS.
        Returns: (count_deleted, log_message)
        """
        if not self.backup_dir.exists():
            return 0, "No backups to clean up"

        now = datetime.utcnow()
        deleted_count = 0
        failed_dirs = []

        for backup_subdir in sorted(self.backup_dir.iterdir()):
            if not backup_subdir.is_dir():
                continue

            # Parse timestamp from directory name (YYYY-MM-DDTHH-MM-SSZ)
            dir_name = backup_subdir.name
            try:
                # Handle both underscore and hyphen-based timestamps
                if "T" in dir_name:
                    # Format: 2026-05-04T14-30-45Z
                    timestamp_str = dir_name.replace("Z", "").replace("-", ":", 2).replace("-", "-", 2)
                    backup_time = datetime.fromisoformat(timestamp_str.replace("T", "T"))
                else:
                    # Format: 2026-05-04 (old format)
                    backup_time = datetime.strptime(dir_name, "%Y-%m-%d")
            except (ValueError, AttributeError):
                # Skip directories with invalid timestamps
                continue

            age = now - backup_time
            if age > timedelta(days=self.BACKUP_RETENTION_DAYS):
                try:
                    shutil.rmtree(backup_subdir)
                    deleted_count += 1
                    self.log_message(f"🗑️  Deleted old backup: {dir_name}")
                except (OSError, PermissionError) as e:
                    failed_dirs.append((dir_name, str(e)))
                    self.log_message(f"⚠️  Warning: Could not clean up {backup_subdir.name} ({e})")

        if deleted_count > 0:
            msg = f"🗑️  Cleaned up {deleted_count} backup(s) older than {self.BACKUP_RETENTION_DAYS} days"
        else:
            msg = "No backups to clean up"

        return deleted_count, msg

    def get_backup_age_warnings(self) -> list:
        """
        Check for backups approaching the 30-day expiration (within 2 days).
        Returns: List of warning messages
        """
        if not self.backup_dir.exists():
            return []

        warnings = []
        now = datetime.utcnow()

        for backup_subdir in sorted(self.backup_dir.iterdir()):
            if not backup_subdir.is_dir():
                continue

            dir_name = backup_subdir.name
            try:
                if "T" in dir_name:
                    timestamp_str = dir_name.replace("Z", "").replace("-", ":", 2).replace("-", "-", 2)
                    backup_time = datetime.fromisoformat(timestamp_str.replace("T", "T"))
                else:
                    backup_time = datetime.strptime(dir_name, "%Y-%m-%d")
            except (ValueError, AttributeError):
                continue

            age = now - backup_time
            days_remaining = self.BACKUP_RETENTION_DAYS - age.days

            # Warn if within 2 days of expiration
            if 0 <= days_remaining <= 2:
                expiration_date = (backup_time + timedelta(days=self.BACKUP_RETENTION_DAYS)).strftime("%Y-%m-%d")
                warnings.append(
                    f"⚠️  Backup from {dir_name} is {age.days} days old "
                    f"(expires in {days_remaining} day{'s' if days_remaining != 1 else ''})"
                )

        return warnings

    # ==== Phase: Audit Trail Logging ====

    def log_audit_trail(
        self,
        status: str,
        commit_sha: str,
        gpg_signature: str,
        backups_cleaned: int,
        files_processed: Optional[int] = None,
        error_msg: Optional[str] = None,
    ) -> str:
        """
        Write audit trail entry to JSON log file.
        Returns: path to log file
        """
        self.log_dir.mkdir(parents=True, exist_ok=True)
        log_path = self.log_dir / self.LOG_FILE

        entry = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "command": "/hyper-update",
            "upstream_repo": "https://github.com/tjmustard/Hypergraph-Coding-Agent-Framework.git",
            "upstream_commit_sha": commit_sha,
            "gpg_signature": gpg_signature,
            "backups_cleaned": backups_cleaned,
            "status": status,
        }

        if files_processed is not None:
            entry["files_processed"] = files_processed

        if error_msg:
            entry["error"] = error_msg

        # Append to log file
        try:
            with open(log_path, "a") as f:
                json.dump(entry, f)
                f.write("\n")

            # Check for log rotation
            self._rotate_logs(log_path)

            return str(log_path)
        except (OSError, IOError) as e:
            self.log_message(f"⚠️  Warning: Could not write to log file ({e})")
            return str(log_path)

    def _rotate_logs(self, log_path: Path):
        """Rotate log file if it exceeds max size; keep last LOG_MAX_FILES files."""
        try:
            if not log_path.exists():
                return

            if log_path.stat().st_size > self.LOG_MAX_SIZE:
                # Rotate files: log.1 -> log.2, etc.
                for i in range(self.LOG_MAX_FILES - 1, 0, -1):
                    old = log_path.parent / f"{log_path.name}.{i}"
                    new = log_path.parent / f"{log_path.name}.{i + 1}"
                    if old.exists():
                        old.replace(new)

                # Move current to .1
                log_path.replace(log_path.parent / f"{log_path.name}.1")

                # Delete old files beyond LOG_MAX_FILES
                for i in range(self.LOG_MAX_FILES + 1, 100):
                    old = log_path.parent / f"{log_path.name}.{i}"
                    if old.exists():
                        old.unlink()
                    else:
                        break

                self.log_message(f"📋 Log rotated (kept last {self.LOG_MAX_FILES} files)")
        except (OSError, PermissionError) as e:
            self.log_message(f"⚠️  Warning: Could not rotate logs ({e})")

    # ==== Utility: Format backup summary ====

    def format_backup_summary(self) -> str:
        """Format a summary of available backups with their ages."""
        if not self.backup_dir.exists() or not list(self.backup_dir.iterdir()):
            return ""

        now = datetime.utcnow()
        summary_lines = ["Available backups:"]

        for backup_subdir in sorted(self.backup_dir.iterdir(), reverse=True)[:10]:  # Last 10
            if not backup_subdir.is_dir():
                continue

            dir_name = backup_subdir.name
            try:
                if "T" in dir_name:
                    timestamp_str = dir_name.replace("Z", "").replace("-", ":", 2).replace("-", "-", 2)
                    backup_time = datetime.fromisoformat(timestamp_str.replace("T", "T"))
                else:
                    backup_time = datetime.strptime(dir_name, "%Y-%m-%d")
            except (ValueError, AttributeError):
                continue

            age = now - backup_time
            days_remaining = self.BACKUP_RETENTION_DAYS - age.days

            age_str = f"{age.days} days old"
            if days_remaining <= 0:
                age_str += " (EXPIRED)"
            elif days_remaining <= 2:
                age_str += f" (expires in {days_remaining} day{'s' if days_remaining != 1 else ''})"

            summary_lines.append(f"- {dir_name} ({age_str})")

        return "\n".join(summary_lines)
