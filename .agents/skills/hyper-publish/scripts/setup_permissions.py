#!/usr/bin/env python3
"""
Patch .claude/settings.json to pre-approve git and gh commands.

Adds or updates the "permissions" key with "allow" patterns for Bash(git *) and Bash(gh *).
This prevents Claude Code from prompting for approval on every git/gh invocation.

Idempotent — safe to run multiple times. Skips if already configured.

Exit code 0 on success (added, updated, or already present).
Exit code 1 on error (e.g., cannot read/write file).
"""

import json
import sys
from pathlib import Path


def patch_settings_json():
    """Patch .claude/settings.json with git/gh permissions."""

    settings_path = Path(".claude/settings.json")

    # Read existing settings
    if settings_path.exists():
        try:
            with open(settings_path, "r") as f:
                settings = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"❌ Error reading .claude/settings.json: {e}", file=sys.stderr)
            return False
    else:
        settings = {}

    # Ensure "permissions" key exists
    if "permissions" not in settings:
        settings["permissions"] = {"allow": []}
    elif not isinstance(settings.get("permissions"), dict):
        print(
            "❌ settings.json 'permissions' is not a dict. Please fix manually.",
            file=sys.stderr,
        )
        return False

    if "allow" not in settings["permissions"]:
        settings["permissions"]["allow"] = []

    # Add git and gh patterns if not already present
    allow_list = settings["permissions"]["allow"]
    git_pattern = "Bash(git *)"
    gh_pattern = "Bash(gh *)"
    touch_pattern = "Bash(touch *)"

    added = []
    for pattern in [git_pattern, gh_pattern, touch_pattern]:
        if pattern not in allow_list:
            allow_list.append(pattern)
            added.append(pattern)

    # Write back if changes were made
    if added:
        try:
            with open(settings_path, "w") as f:
                json.dump(settings, f, indent=2)
                f.write("\n")  # Add trailing newline for good style
            print(f"✅ Permissions updated. Added: {', '.join(added)}", file=sys.stderr)
            return True
        except IOError as e:
            print(
                f"❌ Error writing .claude/settings.json: {e}", file=sys.stderr
            )
            return False
    else:
        print(
            "✅ git/gh permissions already configured in .claude/settings.json",
            file=sys.stderr,
        )
        return True


def main():
    """Main entry point."""
    if patch_settings_json():
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
