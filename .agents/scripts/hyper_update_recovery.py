#!/usr/bin/env python3
"""
hyper_update_recovery.py — Backup management and interrupt recovery for hyper-update

Implements:
- Phase A: Timestamped backup directory creation (ISO 8601)
- Phase B: Backup before mutation (copy files, verify)
- Phase C: Merge state tracking (lock file with JSON state)
- Phase D: Interrupt detection & recovery (resume from last approved section)
- Phase E: /hyper-recover command (list backups, restore files)
"""

import sys
import json
import shutil
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple


class RecoveryManager:
    """Handles backup management, state tracking, and interrupt recovery."""

    BACKUP_BASE_DIR = ".agents/.backup"
    LOCK_FILE = ".agents/.merge-in-progress.lock"
    STALE_LOCK_HOURS = 24

    def __init__(self, project_root: str = "."):
        self.project_root = Path(project_root).resolve()
        self.backup_dir = None  # Created lazily
        self.lock_file = self.project_root / self.LOCK_FILE
        self.current_timestamp = self._get_timestamp()

    def log_message(self, msg: str):
        """Print a message to stderr."""
        print(f"[hyper-update] {msg}", file=sys.stderr)

    def error(self, msg: str):
        """Log an error and exit."""
        self.log_message(f"❌ Error: {msg}")
        sys.exit(1)

    # ==== Phase A: Timestamped Backup Directory Creation ====

    def _get_timestamp(self) -> str:
        """Generate ISO 8601 timestamp in UTC: YYYY-MM-DDTHH-MM-SSZ"""
        return datetime.utcnow().strftime("%Y-%m-%dT%H-%M-%SZ")

    def _create_backup_dir(self) -> Path:
        """
        Create backup directory with current timestamp (lazy creation).
        Only called when first file is actually backed up.
        """
        if self.backup_dir is None:
            backup_path = self.project_root / self.BACKUP_BASE_DIR / self.current_timestamp
            try:
                backup_path.mkdir(parents=True, exist_ok=True)
                self.backup_dir = backup_path
                self.log_message(f"🔐 Backup created: {self.BACKUP_BASE_DIR}/{self.current_timestamp}/")
            except OSError as e:
                self.error(f"Failed to create backup directory: {e}")

        return self.backup_dir

    # ==== Phase B: Backup Before Mutation ====

    def backup_file(self, file_path: str) -> bool:
        """
        Back up a file before mutation.
        Returns: True if successful, False otherwise
        """
        file_path = Path(file_path)
        if not file_path.exists():
            self.error(f"File to backup does not exist: {file_path}")

        backup_dir = self._create_backup_dir()

        # Preserve directory structure
        relative_path = file_path.relative_to(self.project_root)
        backup_path = backup_dir / relative_path

        try:
            backup_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(file_path, backup_path)

            # Verify backup
            if not backup_path.exists() or backup_path.stat().st_size == 0:
                raise OSError(f"Backup verification failed: file size is 0")

            self.log_message(f"Backed up {file_path.name}")
            return True
        except (OSError, IOError) as e:
            self.error(f"Backup failed for {file_path}: {e}")

    # ==== Phase C: Merge State Tracking ====

    def save_merge_state(
        self,
        filename: str,
        sections_total: int,
        sections_approved: int,
        approved_sections: List[str],
        pending_sections: List[str],
    ):
        """Save merge state to lock file for interrupt recovery."""
        lock_data = {
            "start_time": self.current_timestamp,
            "file": filename,
            "sections_total": sections_total,
            "sections_approved": sections_approved,
            "approved_sections": approved_sections,
            "pending_sections": pending_sections,
        }

        try:
            self.lock_file.parent.mkdir(parents=True, exist_ok=True)
            with open(self.lock_file, "w") as f:
                json.dump(lock_data, f, indent=2)
        except (OSError, IOError) as e:
            self.log_message(f"⚠️  Warning: Could not save merge state: {e}")

    def load_merge_state(self) -> Optional[Dict]:
        """Load merge state from lock file."""
        if not self.lock_file.exists():
            return None

        try:
            with open(self.lock_file, "r") as f:
                return json.load(f)
        except (OSError, IOError, json.JSONDecodeError):
            return None

    def clear_merge_state(self):
        """Delete the lock file when merge completes."""
        if self.lock_file.exists():
            try:
                self.lock_file.unlink()
            except OSError:
                self.log_message(f"⚠️  Warning: Could not delete lock file")

    # ==== Phase D: Interrupt Detection & Recovery ====

    def check_for_interrupt(self) -> Optional[Dict]:
        """
        Check if merge was interrupted and prompt user.
        Returns: merge state if user chooses to resume, None if starting fresh
        """
        state = self.load_merge_state()
        if state is None:
            return None

        start_time = state.get("start_time", "unknown")
        lock_age = self._get_lock_age(start_time)

        # Check if lock is stale (older than 24 hours)
        if lock_age > self.STALE_LOCK_HOURS:
            self.log_message(
                f"\n⚠️  Merge in progress from {start_time} (older than 24 hours)"
            )
            response = input("This lock appears stale. Start fresh upgrade? (Y/n): ").lower().strip()
            if response in ("y", "yes", ""):
                self.clear_merge_state()
                return None
            # Otherwise user chose to resume despite staleness

        # Lock is recent; offer to resume
        self.log_message(f"\n⚠️  Merge in progress from {start_time}")
        response = input("Resume merge? (Y/n): ").lower().strip()

        if response in ("y", "yes", ""):
            return state
        else:
            self.clear_merge_state()
            return None

    def _get_lock_age(self, timestamp_str: str) -> int:
        """Calculate age of lock in hours."""
        try:
            lock_time = datetime.fromisoformat(timestamp_str.replace("Z", "+00:00"))
            age = datetime.now(lock_time.tzinfo) - lock_time
            return int(age.total_seconds() / 3600)
        except (ValueError, AttributeError):
            return 0

    # ==== Phase E: /hyper-recover Command ====

    def list_backups(self) -> str:
        """List available backups with dates."""
        backup_base = self.project_root / self.BACKUP_BASE_DIR
        if not backup_base.exists():
            return "No backups available"

        lines = ["Available backups:"]
        backups = sorted(backup_base.iterdir(), reverse=True)

        for backup_dir in backups:
            if not backup_dir.is_dir():
                continue

            dir_name = backup_dir.name
            age_str = self._format_backup_age(dir_name)
            lines.append(f"- {self.BACKUP_BASE_DIR}/{dir_name}/ ({age_str})")

        return "\n".join(lines)

    def _format_backup_age(self, timestamp_str: str) -> str:
        """Format backup age as human-readable string."""
        try:
            if "T" in timestamp_str:
                backup_time = datetime.fromisoformat(timestamp_str.replace("Z", "+00:00"))
            else:
                backup_time = datetime.strptime(timestamp_str, "%Y-%m-%d")

            age = datetime.now(backup_time.tzinfo) - backup_time
            days = age.days

            if days == 0:
                return "today"
            elif days == 1:
                return "1 day old"
            else:
                return f"{days} days old"
        except (ValueError, AttributeError):
            return "unknown age"

    def restore_file(self, filename: str, backup_date: str) -> bool:
        """
        Restore a file from a backup.
        Creates meta-backup of current version before restoring.
        Returns: True if successful, False otherwise
        """
        # Find backup directory matching the date
        backup_base = self.project_root / self.BACKUP_BASE_DIR
        backup_source = None

        for backup_dir in backup_base.iterdir():
            if not backup_dir.is_dir():
                continue
            # Match either full timestamp or just the date part
            if backup_dir.name.startswith(backup_date):
                backup_source = backup_dir / filename
                break

        if backup_source is None or not backup_source.exists():
            self.error(f"Backup not found: {filename} from {backup_date}")

        file_path = self.project_root / filename

        # Create meta-backup of current version
        if file_path.exists():
            meta_backup_timestamp = self._get_timestamp()
            meta_backup_dir = self.project_root / self.BACKUP_BASE_DIR / meta_backup_timestamp
            meta_backup_dir.mkdir(parents=True, exist_ok=True)
            meta_backup_path = meta_backup_dir / filename

            try:
                shutil.copy2(file_path, meta_backup_path)
                self.log_message(f"💾 Meta-backup created: {self.BACKUP_BASE_DIR}/{meta_backup_timestamp}/{filename}")
            except (OSError, IOError) as e:
                self.error(f"Failed to create meta-backup: {e}")

        # Restore from backup
        try:
            file_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(backup_source, file_path)
            self.log_message(
                f"✅ Restored {filename} from {backup_source.parent.name}"
            )
            return True
        except (OSError, IOError) as e:
            self.error(f"Restore failed: {e}")

    def get_available_dates(self) -> List[str]:
        """Get list of available backup dates."""
        backup_base = self.project_root / self.BACKUP_BASE_DIR
        if not backup_base.exists():
            return []

        dates = set()
        for backup_dir in backup_base.iterdir():
            if backup_dir.is_dir():
                # Extract just the date part (YYYY-MM-DD)
                timestamp = backup_dir.name
                if "T" in timestamp:
                    date_part = timestamp.split("T")[0]
                    dates.add(date_part)
                else:
                    dates.add(timestamp)

        return sorted(dates, reverse=True)
