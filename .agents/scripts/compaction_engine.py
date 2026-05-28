#!/usr/bin/env python3
"""
Context Compaction Engine

Collapses full /hyper-resolve debate transcripts into "Final Truth" documents
that preserve trade-off reasoning while discarding verbatim debate text.
"""

import hashlib
import re
from datetime import datetime
from dataclasses import dataclass
from typing import List, Optional, Dict, Any


@dataclass
class TradeoffEntry:
    """Represents a single decision made during resolution."""
    title: str
    architect_position: str
    red_team_position: str
    resolution: str
    risks: List[str]
    mitigations: List[str]


@dataclass
class CompactionResult:
    """Result of compaction operation."""
    success: bool
    final_truth_content: str
    archived_transcript_path: str
    tradeoff_count: int
    error_message: Optional[str] = None


class CompactionEngine:
    """Extracts trade-off decisions from /hyper-resolve transcripts."""

    def __init__(self, feature_name: str = "feature"):
        self.feature_name = feature_name
        self.tradeoffs: List[TradeoffEntry] = []

    def compact(self, transcript: str) -> CompactionResult:
        """
        Extract trade-off entries from a debate transcript.

        Args:
            transcript: Full /hyper-resolve debate text

        Returns:
            CompactionResult with final_truth content and archive path
        """
        try:
            # Parse tradeoffs from transcript
            self.tradeoffs = self._extract_tradeoffs(transcript)

            if not self.tradeoffs:
                return CompactionResult(
                    success=False,
                    final_truth_content="",
                    archived_transcript_path="",
                    tradeoff_count=0,
                    error_message="No tradeoffs found in transcript"
                )

            # Generate archive path
            transcript_hash = self._hash_transcript(transcript)
            archive_path = self._generate_archive_path(transcript_hash)

            # Generate Final Truth document
            final_truth = self._generate_final_truth(archive_path)

            return CompactionResult(
                success=True,
                final_truth_content=final_truth,
                archived_transcript_path=archive_path,
                tradeoff_count=len(self.tradeoffs),
                error_message=None
            )
        except Exception as e:
            return CompactionResult(
                success=False,
                final_truth_content="",
                archived_transcript_path="",
                tradeoff_count=0,
                error_message=f"Compaction failed: {str(e)}"
            )

    def _extract_tradeoffs(self, transcript: str) -> List[TradeoffEntry]:
        """
        Parse transcript to extract tradeoff entries.
        Looks for structured sections with Architect/Red Team/Resolution.
        """
        tradeoffs = []

        # Pattern: match decision blocks with three clear positions
        # This is a simplified regex pattern; real implementation would be more robust
        decision_pattern = r"(?:## |###\s+)(?:Decision|Tradeoff):\s*(.+?)\n\n(.*?)(?=(?:## |###\s+(?:Decision|Tradeoff)|$))"

        for match in re.finditer(decision_pattern, transcript, re.DOTALL | re.IGNORECASE):
            title = match.group(1).strip()
            content = match.group(2)

            architect_pos = self._extract_section(content, "architect")
            red_team_pos = self._extract_section(content, "red team|red_team")
            resolution = self._extract_section(content, "resolution|final decision")
            risks = self._extract_list_items(content, "risk")
            mitigations = self._extract_list_items(content, "mitigation")

            if architect_pos and red_team_pos and resolution:
                tradeoffs.append(TradeoffEntry(
                    title=title,
                    architect_position=architect_pos,
                    red_team_position=red_team_pos,
                    resolution=resolution,
                    risks=risks,
                    mitigations=mitigations
                ))

        return tradeoffs

    def _extract_section(self, text: str, section_name: str) -> str:
        """Extract content from a named section."""
        pattern = rf"(?:{section_name}):\s*(.+?)(?=\n\n|\n[A-Z]|$)"
        match = re.search(pattern, text, re.IGNORECASE | re.DOTALL)
        return match.group(1).strip() if match else ""

    def _extract_list_items(self, text: str, list_name: str) -> List[str]:
        """Extract list items from a section (e.g., risks, mitigations)."""
        pattern = rf"(?:{list_name}s?):\s*(.*?)(?=\n\n|\n[A-Z]|$)"
        match = re.search(pattern, text, re.IGNORECASE | re.DOTALL)
        if not match:
            return []

        items_text = match.group(1)
        items = [line.strip() for line in items_text.split("\n") if line.strip().startswith(("-", "*", "•"))]
        return [item.lstrip("-*•").strip() for item in items]

    def _hash_transcript(self, transcript: str) -> str:
        """Generate 8-char hash of transcript for collision prevention."""
        sha256 = hashlib.sha256(transcript.encode())
        return sha256.hexdigest()[:8]

    def _generate_archive_path(self, transcript_hash: str) -> str:
        """Generate archive file path with date and hash."""
        date_str = datetime.utcnow().strftime("%Y-%m-%d")
        return f"spec/archive/{date_str}_{transcript_hash}_debate.md"

    def _generate_final_truth(self, archive_path: str) -> str:
        """Generate the Final Truth markdown document."""
        timestamp = datetime.utcnow().isoformat()

        doc = f"""# Final Truth: {self.feature_name}

**Generated:** {timestamp}
**Feature:** {self.feature_name}
**Archive:** See full transcript: `{archive_path}`

---

## Trade-off Decisions

"""

        for i, tradeoff in enumerate(self.tradeoffs, 1):
            doc += f"""### {i}. {tradeoff.title}

**Architect Position:**
{tradeoff.architect_position}

**Red Team Position:**
{tradeoff.red_team_position}

**Resolution:**
{tradeoff.resolution}

"""
            if tradeoff.risks:
                doc += "**Risks Identified:**\n"
                for risk in tradeoff.risks:
                    doc += f"- {risk}\n"
                doc += "\n"

            if tradeoff.mitigations:
                doc += "**Mitigations:**\n"
                for mitigation in tradeoff.mitigations:
                    doc += f"- {mitigation}\n"
                doc += "\n"

        doc += f"""---

## Archive & Retention

Full debate transcript archived at: `{archive_path}`

**Retention Policy:** Minimum 90 days. This transcript is the authoritative record
of all decisions made and rejected during this feature cycle.

**Link back to compiled specs:** See `spec/compiled/` for SuperPRD and MiniPRDs.
"""

        return doc


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python compaction_engine.py <path_to_transcript> [feature_name]")
        sys.exit(1)

    transcript_path = sys.argv[1]
    feature_name = sys.argv[2] if len(sys.argv) > 2 else "feature"

    try:
        with open(transcript_path, 'r') as f:
            transcript = f.read()

        engine = CompactionEngine(feature_name)
        result = engine.compact(transcript)

        if result.success:
            print(f"✅ COMPACTION SUCCESS")
            print(f"Tradeoffs extracted: {result.tradeoff_count}")
            print(f"Archive path: {result.archived_transcript_path}")
            print(f"\nFinal Truth (first 500 chars):\n{result.final_truth_content[:500]}...\n")
        else:
            print(f"❌ COMPACTION FAILED: {result.error_message}")
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Transcript file not found: {transcript_path}")
        sys.exit(1)
