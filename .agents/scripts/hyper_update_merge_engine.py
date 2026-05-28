#!/usr/bin/env python3
"""
hyper_update_merge_engine.py — Sensitive file merge orchestration with per-section approval

Implements:
- File interaction prompting (Replace/Merge/Skip)
- Markdown section parsing by ### headings (preserves code blocks)
- Per-section approval loop with diff display
- File assembly and validation (YAML, markdown)
- Summary report generation
"""

import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import yaml


@dataclass
class Section:
    """Represents a markdown section."""
    heading: str  # "### Section Name" or "" for preamble
    content: str  # Content of this section
    line_num: int  # Starting line number


class MergeEngine:
    """Handles section-by-section merging of sensitive files."""

    def __init__(self, project_root: str = "."):
        self.project_root = Path(project_root).resolve()
        self.changes = {
            "auto_updated": [],
            "replaced": [],
            "merged": [],
            "skipped": [],
        }

    def log_message(self, msg: str):
        """Print a message to stderr for real-time feedback."""
        print(f"[hyper-update] {msg}", file=sys.stderr)

    def error(self, msg: str):
        """Log an error and exit."""
        self.log_message(f"❌ Error: {msg}")
        sys.exit(1)

    # ==== Phase A: File Interaction Prompting ====

    def prompt_file_action(self, filename: str, diff_preview: str) -> str:
        """
        Prompt user for action on a file with changes.
        Returns: "replace", "merge", or "skip"
        """
        self.log_message(f"\n📝 File: {filename}")
        self.log_message("Diff preview (first 10 lines):")
        for line in diff_preview.split("\n")[:10]:
            self.log_message(f"  {line}")

        while True:
            response = input(
                "\n[R]eplace / [M]erge / [S]kip? (r/m/s): "
            ).lower().strip()

            if response in ("r", "replace"):
                return "replace"
            elif response in ("m", "merge"):
                return "merge"
            elif response in ("s", "skip"):
                return "skip"
            else:
                self.log_message("Invalid input. Enter r, m, or s.")

    def generate_diff_preview(self, upstream: str, local: str) -> str:
        """Generate a unified diff-style preview (simplified)."""
        upstream_lines = upstream.split("\n")
        local_lines = local.split("\n")

        diff_lines = []
        for i, (up, loc) in enumerate(zip(upstream_lines[:5], local_lines[:5])):
            if up != loc:
                diff_lines.append(f"- {up[:60]}")
                diff_lines.append(f"+ {loc[:60]}")
            else:
                diff_lines.append(f"  {up[:60]}")

        return "\n".join(diff_lines)

    # ==== Phase B: Markdown Section Parsing ====

    def parse_sections(self, content: str) -> List[Section]:
        """
        Parse markdown into sections by ### headings.
        Preserves code blocks (does not split within triple backticks).
        Returns: List of Section objects
        """
        lines = content.split("\n")
        sections = []
        current_section_heading = ""
        current_section_lines = []
        in_code_block = False

        for i, line in enumerate(lines):
            # Track code block state
            if line.strip().startswith("```"):
                in_code_block = not in_code_block

            # Check for section heading (only outside code blocks)
            if line.startswith("###") and not in_code_block:
                # Save previous section
                if current_section_lines or current_section_heading:
                    sections.append(
                        Section(
                            heading=current_section_heading,
                            content="\n".join(current_section_lines),
                            line_num=len(sections),
                        )
                    )

                current_section_heading = line
                current_section_lines = []
            else:
                current_section_lines.append(line)

        # Save final section
        if current_section_lines or current_section_heading:
            sections.append(
                Section(
                    heading=current_section_heading,
                    content="\n".join(current_section_lines),
                    line_num=len(sections),
                )
            )

        return sections

    # ==== Phase C: Per-Section Approval Loop ====

    def approve_sections(
        self, upstream_sections: List[Section], local_sections: List[Section]
    ) -> Dict[str, str]:
        """
        Show each differing section and prompt for approval.
        Returns: Dict of heading -> approved content
        """
        approved = {}
        local_by_heading = {s.heading: s.content for s in local_sections}
        upstream_by_heading = {s.heading: s.content for s in upstream_sections}

        # Process each upstream section
        total_sections = len(upstream_sections)
        approved_count = 0

        for i, upstream_sec in enumerate(upstream_sections):
            heading = upstream_sec.heading
            local_content = local_by_heading.get(heading, None)

            # If section unchanged, approve silently
            if local_content is not None and local_content == upstream_sec.content:
                approved[heading] = local_content
                approved_count += 1
                continue

            # Section differs or is new; prompt user
            self.log_message(f"\n📋 Section {approved_count + 1} of {total_sections}: {heading}")

            if local_content is None:
                self.log_message("ℹ️  New section in upstream (not in local)")
                self.log_message(f"Content:\n{upstream_sec.content[:500]}...")
                response = input("Accept upstream section? (y/n): ").lower().strip()
                if response in ("y", "yes"):
                    approved[heading] = upstream_sec.content
                else:
                    # Skip new section
                    pass
                approved_count += 1
            else:
                # Section exists in both; show diff
                self.log_message("Upstream version:")
                self.log_message(upstream_sec.content[:300])
                self.log_message("\nLocal version:")
                self.log_message(local_content[:300])

                while True:
                    response = input(
                        "\n[K]eep local / [A]ccept upstream / [E]dit manually? (k/a/e): "
                    ).lower().strip()

                    if response in ("k", "keep"):
                        approved[heading] = local_content
                        break
                    elif response in ("a", "accept"):
                        approved[heading] = upstream_sec.content
                        break
                    elif response in ("e", "edit"):
                        edited = input("Paste edited content:\n")
                        approved[heading] = edited
                        break
                    else:
                        self.log_message("Invalid input. Enter k, a, or e.")

                approved_count += 1

        return approved

    # ==== Phase D: File Assembly & Validation ====

    def assemble_file(self, approved_sections: Dict[str, str]) -> str:
        """Assemble approved sections into final file content."""
        parts = []
        for heading, content in approved_sections.items():
            if heading:  # Include heading if present
                parts.append(heading)
            parts.append(content)

        return "\n".join(parts)

    def validate_file(self, filename: str, content: str) -> Tuple[bool, str]:
        """
        Validate merged file for syntax errors.
        Returns: (is_valid, error_message)
        """
        if filename.endswith(".yml") or filename.endswith(".yaml"):
            try:
                yaml.safe_load(content)
                return True, ""
            except yaml.YAMLError as e:
                return False, f"YAML syntax error: {e}"

        if filename.endswith(".md"):
            # Check for unmatched code block delimiters
            backtick_count = content.count("```")
            if backtick_count % 2 != 0:
                return (
                    False,
                    "Markdown has unmatched code block delimiters (```)",
                )

        return True, ""

    def handle_validation_failure(self, filename: str, error: str) -> str:
        """Prompt user after validation failure."""
        self.log_message(f"\n❌ Validation failed: {error}")

        while True:
            response = input(
                "What would you like to do?\n[K]eep original / [R]e-merge / [A]ccept anyway? (k/r/a): "
            ).lower().strip()

            if response in ("k", "keep"):
                return "keep_original"
            elif response in ("r", "remerge"):
                return "remerge"
            elif response in ("a", "accept"):
                return "accept_anyway"
            else:
                self.log_message("Invalid input. Enter k, r, or a.")

    # ==== Phase E: Summary Report ====

    def generate_summary(self, backup_path: Optional[str] = None) -> str:
        """Generate upgrade summary report."""
        lines = ["\n📊 Upgrade Summary"]
        lines.append(f"✅ Auto-updated: {len(self.changes['auto_updated'])} files")
        if self.changes["auto_updated"]:
            for f in self.changes["auto_updated"][:3]:
                lines.append(f"   - {f}")
            if len(self.changes["auto_updated"]) > 3:
                lines.append(f"   ... and {len(self.changes['auto_updated']) - 3} more")

        lines.append(f"✅ Replaced: {len(self.changes['replaced'])} file(s)")
        for f in self.changes["replaced"]:
            lines.append(f"   - {f}")

        lines.append(f"✅ Merged: {len(self.changes['merged'])} file(s)")
        for f in self.changes["merged"]:
            lines.append(f"   - {f}")

        lines.append(f"⏭️  Skipped: {len(self.changes['skipped'])} file(s)")
        for f in self.changes["skipped"]:
            lines.append(f"   - {f}")

        if backup_path:
            lines.append(f"\n🔐 Backup: {backup_path}")

        return "\n".join(lines)

    # ==== Main Orchestration ====

    def merge_file(
        self,
        filename: str,
        upstream_content: str,
        local_content: str,
        backup_dir: Optional[Path] = None,
    ) -> Tuple[bool, Optional[str]]:
        """
        Orchestrate merge for a single file.
        Returns: (success, final_content or None)
        """
        # Phase A: Determine action
        action = self.prompt_file_action(
            filename, self.generate_diff_preview(upstream_content, local_content)
        )

        if action == "skip":
            self.changes["skipped"].append(filename)
            return True, None

        if action == "replace":
            if backup_dir:
                self._backup_file(filename, local_content, backup_dir)
            self.changes["replaced"].append(filename)
            return True, upstream_content

        # action == "merge"
        # Phase B: Parse sections
        upstream_sections = self.parse_sections(upstream_content)
        local_sections = self.parse_sections(local_content)

        # Phase C: Approve sections
        approved = self.approve_sections(upstream_sections, local_sections)

        # Phase D: Assemble and validate
        final_content = self.assemble_file(approved)
        is_valid, error_msg = self.validate_file(filename, final_content)

        if not is_valid:
            action = self.handle_validation_failure(filename, error_msg)
            if action == "keep_original":
                self.changes["skipped"].append(filename)
                return True, None
            elif action == "remerge":
                # Recursively re-merge
                return self.merge_file(
                    filename, upstream_content, local_content, backup_dir
                )
            elif action == "accept_anyway":
                pass  # Proceed with writing

        if backup_dir:
            self._backup_file(filename, local_content, backup_dir)

        self.changes["merged"].append(filename)
        self.log_message(f"✅ Merged {filename} (validation: pass)")
        return True, final_content

    def _backup_file(self, filename: str, content: str, backup_dir: Path):
        """Create backup of original file."""
        backup_path = backup_dir / filename
        backup_path.parent.mkdir(parents=True, exist_ok=True)
        with open(backup_path, "w") as f:
            f.write(content)
