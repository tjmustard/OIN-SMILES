#!/usr/bin/env python3
"""
Override Manager

Manages skill model assignment overrides with scope lifecycle (single_run vs permanent).
Handles override creation, persistence, and cleanup.
"""

import os
import yaml
from datetime import datetime
from typing import Optional


class OverrideManager:
    """Manages model and ceiling overrides for skills."""

    def __init__(self, skills_base_dir: str = ".agents/skills"):
        self.skills_base_dir = skills_base_dir

    def set_override(
        self,
        skill_name: str,
        model: str,
        scope: str = "single_run",
        reason: str = "",
        verbose: bool = True
    ) -> bool:
        """
        Set a model override for a skill.

        Args:
            skill_name: Name of skill (e.g., "hyper-execute")
            model: Target model tier ("opus", "sonnet", "haiku")
            scope: "single_run" or "permanent"
            reason: Reason for override (logged)
            verbose: Print warnings/messages

        Returns:
            Success status
        """
        if model not in ["opus", "sonnet", "haiku"]:
            if verbose:
                print(f"❌ Invalid model: {model}. Must be: opus, sonnet, haiku")
            return False

        if scope not in ["single_run", "permanent"]:
            if verbose:
                print(f"❌ Invalid scope: {scope}. Must be: single_run, permanent")
            return False

        meta_path = self._resolve_meta_path(skill_name)

        # Load existing metadata
        try:
            with open(meta_path, 'r') as f:
                data = yaml.safe_load(f) or {}
        except FileNotFoundError:
            if verbose:
                print(f"⚠️  META.yml not found for {skill_name}. Creating new file.")
            data = {"name": skill_name, "assigned_model": "opus"}
        except Exception as e:
            if verbose:
                print(f"❌ Failed to read META.yml: {e}")
            return False

        # Add override metadata
        data["override"] = {
            "model": model,
            "scope": scope,
            "set_by": "override_manager",
            "set_at": datetime.utcnow().isoformat(),
            "reason": reason
        }

        # Save updated metadata
        try:
            os.makedirs(os.path.dirname(meta_path), exist_ok=True)
            with open(meta_path, 'w') as f:
                yaml.dump(data, f, default_flow_style=False, sort_keys=False)

            if verbose:
                msg = f"⚠️  Override set: {skill_name} → {model} (scope: {scope})"
                if reason:
                    msg += f". Reason: {reason}"
                print(msg)

            return True

        except Exception as e:
            if verbose:
                print(f"❌ Failed to save override: {e}")
            return False

    def clear_override(self, skill_name: str, verbose: bool = True) -> bool:
        """
        Clear an active override for a skill.

        Args:
            skill_name: Name of skill
            verbose: Print messages

        Returns:
            Success status
        """
        meta_path = self._resolve_meta_path(skill_name)

        if not os.path.exists(meta_path):
            if verbose:
                print(f"ℹ️  No META.yml found for {skill_name}. No override to clear.")
            return True

        try:
            with open(meta_path, 'r') as f:
                data = yaml.safe_load(f) or {}

            if "override" in data:
                del data["override"]
                with open(meta_path, 'w') as f:
                    yaml.dump(data, f, default_flow_style=False, sort_keys=False)

                if verbose:
                    print(f"✅ Override cleared for {skill_name}")

            return True

        except Exception as e:
            if verbose:
                print(f"❌ Failed to clear override: {e}")
            return False

    def check_single_run_override(self, skill_name: str) -> bool:
        """
        Check if skill has an active single_run override.

        Returns:
            True if single_run override is active, False otherwise
        """
        meta_path = self._resolve_meta_path(skill_name)

        if not os.path.exists(meta_path):
            return False

        try:
            with open(meta_path, 'r') as f:
                data = yaml.safe_load(f) or {}

            override = data.get("override", {})
            return override.get("scope") == "single_run"

        except Exception:
            return False

    def _resolve_meta_path(self, skill_name: str) -> str:
        """Resolve path to skill's META.yml file."""
        candidates = [
            os.path.join(self.skills_base_dir, skill_name, "META.yml"),
            os.path.join(self.skills_base_dir, skill_name.replace("-", "_"), "META.yml"),
            os.path.join(self.skills_base_dir, f"hyper-{skill_name}", "META.yml"),
        ]
        return candidates[0]


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: python override_manager.py <skill_name> <model> [--scope single_run|permanent] [--reason 'text']")
        print("Example: python override_manager.py hyper-execute sonnet --scope single_run")
        sys.exit(1)

    skill_name = sys.argv[1]
    model = sys.argv[2]
    scope = "single_run"
    reason = ""

    if "--scope" in sys.argv:
        idx = sys.argv.index("--scope")
        if idx + 1 < len(sys.argv):
            scope = sys.argv[idx + 1]

    if "--reason" in sys.argv:
        idx = sys.argv.index("--reason")
        if idx + 1 < len(sys.argv):
            reason = sys.argv[idx + 1]

    manager = OverrideManager()
    success = manager.set_override(skill_name, model, scope, reason, verbose=True)

    sys.exit(0 if success else 1)
