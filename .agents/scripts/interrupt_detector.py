#!/usr/bin/env python3
"""
Interrupt Detector & Recovery

Detects if context compaction was interrupted (e.g., /hyper-resolve killed before /compact)
and automatically triggers recovery compaction at /hyper-execute start.
"""

import os
import glob
from datetime import datetime, timedelta
from typing import Tuple, Optional
from compaction_engine import CompactionEngine


class InterruptDetector:
    """Detects and recovers from incomplete compaction."""

    FINAL_TRUTH_GLOB = "spec/compiled/final_truth_*.md"
    TRANSCRIPT_GLOB = "spec/active/*_debate.md"

    def __init__(self, verbose: bool = True):
        self.verbose = verbose

    def detect_and_recover(self) -> Tuple[bool, str]:
        """
        Check if compaction is missing and recover if needed.

        Returns:
            (recovery_needed: bool, status_message: str)
        """
        # Check if Final Truth exists from recent /hyper-resolve
        final_truth_files = glob.glob(self.FINAL_TRUTH_GLOB)

        if final_truth_files:
            # Sort by modification time; use most recent
            final_truth_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            most_recent = final_truth_files[0]

            # Check if it's fresh (within last 2 hours, assuming /hyper-resolve took < 2h)
            mtime = os.path.getmtime(most_recent)
            age_seconds = (datetime.utcnow().timestamp() - mtime)
            is_fresh = age_seconds < (2 * 3600)  # 2 hours

            if is_fresh:
                msg = f"Compaction file exists and is fresh: {most_recent}"
                if self.verbose:
                    print(f"✅ {msg}")
                return False, msg

        # Compaction missing or stale; attempt recovery
        msg = "Compaction file missing or stale. Attempting auto-recovery..."
        if self.verbose:
            print(f"⚠️  {msg}")

        transcript_files = glob.glob(self.TRANSCRIPT_GLOB)
        if not transcript_files:
            error_msg = "No debate transcript found in spec/active/. Cannot auto-recover."
            if self.verbose:
                print(f"❌ {error_msg}")
            return True, error_msg

        # Use most recent transcript
        transcript_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
        transcript_path = transcript_files[0]

        # Run compaction on the transcript
        try:
            with open(transcript_path, 'r') as f:
                transcript = f.read()

            feature_name = self._extract_feature_name(transcript_path)
            engine = CompactionEngine(feature_name)
            result = engine.compact(transcript)

            if result.success:
                recovery_msg = (
                    f"Auto-recovery successful. Extracted {result.tradeoff_count} "
                    f"tradeoff(s). Archive: {result.archived_transcript_path}"
                )
                if self.verbose:
                    print(f"✅ {recovery_msg}")
                return True, recovery_msg
            else:
                error_msg = f"Auto-recovery failed: {result.error_message}"
                if self.verbose:
                    print(f"❌ {error_msg}")
                return True, error_msg

        except Exception as e:
            error_msg = f"Exception during auto-recovery: {str(e)}"
            if self.verbose:
                print(f"❌ {error_msg}")
            return True, error_msg

    def _extract_feature_name(self, file_path: str) -> str:
        """Extract feature name from transcript file path."""
        # Format: spec/active/YYYYMMDD_HHMMSS_FeatureName_debate.md
        basename = os.path.basename(file_path)
        # Remove _debate.md suffix and timestamp prefix
        parts = basename.replace("_debate.md", "").split("_")
        if len(parts) > 2:
            return "_".join(parts[2:])  # Skip date and time
        return "feature"

    def cleanup_old_archives(self, retention_days: int = 90) -> Tuple[int, str]:
        """
        Clean up archived debate transcripts older than retention period.

        Returns:
            (count_deleted: int, status_message: str)
        """
        archive_files = glob.glob("spec/archive/*_debate.md")
        cutoff_time = datetime.utcnow().timestamp() - (retention_days * 86400)
        count_deleted = 0

        for archive_file in archive_files:
            mtime = os.path.getmtime(archive_file)
            if mtime < cutoff_time:
                try:
                    os.remove(archive_file)
                    count_deleted += 1
                    if self.verbose:
                        print(f"Deleted (beyond {retention_days}-day retention): {archive_file}")
                except Exception as e:
                    if self.verbose:
                        print(f"Failed to delete {archive_file}: {str(e)}")

        msg = f"Cleaned {count_deleted} archived debate transcript(s) beyond {retention_days}-day retention."
        return count_deleted, msg


if __name__ == "__main__":
    import sys

    detector = InterruptDetector(verbose=True)

    # Check for recovery needs
    recovery_needed, status = detector.detect_and_recover()
    print(f"\nRecovery Needed: {recovery_needed}")
    print(f"Status: {status}")

    # Optionally cleanup old archives
    if len(sys.argv) > 1 and sys.argv[1] == "--cleanup":
        count, msg = detector.cleanup_old_archives(retention_days=90)
        print(f"\n{msg}")
