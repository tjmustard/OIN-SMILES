#!/usr/bin/env python3
"""
Generate a commit message suggestion from CHANGELOG [Unreleased] and git diff.

Reads:
  - CHANGELOG.md — extracts [Unreleased] block and recent version blocks for style reference
  - git log — infers commit style (Release vX.Y.Z, feat:, imperative, etc.)
  - git diff — infers scope from changed file names

Outputs JSON:
  {"subject": "...", "body": "...", "source_bullets": [...]}

Exit code 0 on success (even if CHANGELOG is empty — will output fallback).
Exit code 1 on error (e.g., no git repo).
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd, cwd="."):
    """Run shell command and return output, or None on error."""
    try:
        result = subprocess.run(
            cmd, shell=True, cwd=cwd, capture_output=True, text=True, timeout=5
        )
        return result.stdout.strip() if result.returncode == 0 else None
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def parse_changelog_section(section_text):
    """Parse a changelog section block into categories.

    Args: section_text (str) — the content between ## headers

    Returns: dict with 'added', 'changed', 'removed' keys, each a list of bullets.
    """
    result = {"added": [], "changed": [], "removed": []}

    current_section = None
    for line in section_text.split("\n"):
        if re.match(r"^### (Added|Changed|Removed)", line):
            current_section = line.split("### ")[1].lower()
        elif line.startswith("- ") and current_section:
            # Remove markdown bold markers and clean up
            bullet = re.sub(r"\*\*([^*]+)\*\*", r"\1", line[2:].strip())
            # Remove trailing period if present
            bullet = bullet.rstrip(".")
            if bullet:
                result[current_section].append(bullet)

    return result


def parse_changelog_best_block(changelog_path):
    """Extract [Unreleased] block, or fall back to most recent versioned block.

    After hyper-document bumps a version, [Unreleased] becomes empty. This function
    prefers [Unreleased] when non-empty, falls back to the most recent version block.

    Returns: tuple (data_dict, version_string_or_None)
      - data_dict: {'added': [...], 'changed': [...], 'removed': [...]}
      - version_string: '0.4.2' if sourced from [0.4.2], None if from [Unreleased]
    """
    if not Path(changelog_path).exists():
        return {"added": [], "changed": [], "removed": []}, None

    try:
        with open(changelog_path, "r") as f:
            content = f.read()

        # Try [Unreleased] first
        unreleased_match = re.search(
            r"## \[Unreleased\](.*?)(?=## \[|$)", content, re.DOTALL
        )
        if unreleased_match:
            unreleased_block = unreleased_match.group(1)
            result = parse_changelog_section(unreleased_block)
            # Check if any bullets were found
            if any(result.values()):
                return result, None

        # [Unreleased] is empty or doesn't exist — try most recent version block
        version_match = re.search(
            r"## \[(\d+\.\d+\.\d+)\](.*?)(?=## \[|$)", content, re.DOTALL
        )
        if version_match:
            version = version_match.group(1)
            version_block = version_match.group(2)
            result = parse_changelog_section(version_block)
            if any(result.values()):
                return result, version

        # No content found anywhere
        return {"added": [], "changed": [], "removed": []}, None

    except Exception:
        return {"added": [], "changed": [], "removed": []}, None


def infer_style_from_recent_commits():
    """Analyze recent commits to infer style (Release, feat:, imperative, etc)."""
    logs = run_cmd("git log --oneline -10")
    if not logs:
        return "imperative"  # Default fallback

    commit_lines = logs.split("\n")

    release_count = sum(1 for line in commit_lines if "Release v" in line)
    feat_count = sum(1 for line in commit_lines if line.split(None, 1)[1].startswith("feat:") if len(line.split(None, 1)) > 1)

    # If recent commits are mostly releases, suggest Release format
    if release_count > 3:
        return "release"
    # If feat: prefix is common, suggest Conventional Commits
    elif feat_count > 3:
        return "conventional"
    # Otherwise, imperative (bare verb)
    return "imperative"


def extract_changed_paths():
    """Get list of changed files for scope hints."""
    paths = run_cmd("git diff --name-only HEAD")
    if not paths:
        return []
    return paths.split("\n")


def generate_subject(changelog_data, style, changed_paths, version=None):
    """Generate a commit subject line based on changelog and style.

    Args:
        version (str): If provided (e.g., "0.4.2"), use Release format regardless of style.
    """

    # Count changes by type
    added = len(changelog_data.get("added", []))
    changed = len(changelog_data.get("changed", []))
    removed = len(changelog_data.get("removed", []))

    total = added + changed + removed

    # Fallback if changelog is empty: infer from file paths
    if total == 0:
        if changed_paths:
            # Group by directory
            prefixes = [p.split("/")[0] for p in changed_paths]
            prefix = max(set(prefixes), key=prefixes.count) if prefixes else "files"
            return f"Update {prefix}"
        return "Update files"

    # Build subject based on dominant change type
    if added > changed and added > removed:
        dominant = "Added"
        first_item = changelog_data["added"][0] if changelog_data["added"] else ""
    elif removed > changed:
        dominant = "Removed"
        first_item = changelog_data["removed"][0] if changelog_data["removed"] else ""
    else:
        dominant = "Changed"
        first_item = changelog_data["changed"][0] if changelog_data["changed"] else ""

    # If version is provided, use Release format
    if version:
        return f"Release v{version}: {first_item[:50]}"

    # Format subject based on style preference
    if style == "release":
        return f"{dominant}: {first_item[:50]}"
    elif style == "conventional":
        action = "feat" if "added" in dominant.lower() else "fix" if "changed" in dominant.lower() else "refactor"
        return f"{action}: {first_item[:50]}"
    else:  # imperative
        verb_map = {"Added": "Add", "Changed": "Update", "Removed": "Remove"}
        verb = verb_map.get(dominant, "Update")
        return f"{verb} {first_item[:50]}"


def generate_body(changelog_data):
    """Generate optional commit body from multiple changelog bullets."""
    # If there are multiple bullets, include them in body
    all_bullets = (
        changelog_data.get("added", []) +
        changelog_data.get("changed", []) +
        changelog_data.get("removed", [])
    )

    if len(all_bullets) <= 1:
        return ""  # Single item; no body needed

    # Format remaining bullets as body
    body_bullets = [f"  - {b}" for b in all_bullets[:5]]  # Limit to 5 bullets
    return "\n".join(body_bullets)


def main():
    """Main entry point."""
    try:
        # Verify git repo
        if not run_cmd("git rev-parse --git-dir"):
            print(json.dumps({"error": "not a git repository"}), file=sys.stderr)
            sys.exit(1)

        # Read CHANGELOG and git state
        changelog_data, version = parse_changelog_best_block("CHANGELOG.md")
        style = infer_style_from_recent_commits()
        changed_paths = extract_changed_paths()

        # Generate subject and body
        subject = generate_subject(changelog_data, style, changed_paths, version=version)
        body = generate_body(changelog_data)

        # Ensure subject is ≤72 chars
        subject = subject[:72]

        # Compile result
        result = {
            "subject": subject,
            "body": body,
            "style": style if not version else "release",
            "source_bullets": (
                changelog_data.get("added", []) +
                changelog_data.get("changed", []) +
                changelog_data.get("removed", [])
            ),
        }

        print(json.dumps(result))
        sys.exit(0)

    except Exception as e:
        print(json.dumps({"error": str(e)}), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
