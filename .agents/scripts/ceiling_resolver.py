#!/usr/bin/env python3
"""
Thinking Token Ceiling Resolver

Determines the appropriate thinking token ceiling for a skill based on
model tier and per-skill overrides.
"""

import os
import yaml
from typing import Optional


class CeilingResolver:
    """Resolves thinking token ceiling for skills."""

    # Default thinking token ceilings by model tier
    MODEL_DEFAULTS = {
        "opus": 20000,
        "sonnet": 10000,
        "haiku": 2000,
    }

    def __init__(self, skills_base_dir: str = ".agents/skills"):
        self.skills_base_dir = skills_base_dir

    def get_ceiling(
        self,
        skill_name: str,
        model: str = "opus",
        verbose: bool = True
    ) -> int:
        """
        Get the thinking token ceiling for a skill.

        Args:
            skill_name: Name of skill (e.g., "hyper-execute")
            model: Model tier ("opus", "sonnet", "haiku")
            verbose: Log ceiling source

        Returns:
            Thinking token ceiling (integer)
        """
        # Try to load skill-specific override
        ceiling = self._load_skill_override(skill_name)

        if ceiling is not None:
            if verbose:
                print(f"ℹ️  {skill_name}: Using custom ceiling ({ceiling} tokens)")
            return ceiling

        # Fall back to model tier default
        default = self.MODEL_DEFAULTS.get(model, 2000)

        if verbose:
            print(f"ℹ️  {skill_name}: Using {model} default ceiling ({default} tokens)")

        return default

    def set_override(
        self,
        skill_name: str,
        tokens: int,
        reason: str = "",
        verbose: bool = True
    ) -> bool:
        """
        Set a custom thinking token ceiling for a skill.

        Args:
            skill_name: Name of skill
            tokens: Token ceiling (must be positive)
            reason: Reason for override
            verbose: Log messages

        Returns:
            Success status
        """
        if tokens <= 0:
            if verbose:
                print(f"❌ Token ceiling must be positive. Got: {tokens}")
            return False

        meta_path = self._resolve_meta_path(skill_name)

        # Load existing metadata
        try:
            with open(meta_path, 'r') as f:
                data = yaml.safe_load(f) or {}
        except FileNotFoundError:
            data = {"name": skill_name, "assigned_model": "opus"}
        except Exception as e:
            if verbose:
                print(f"❌ Failed to read META.yml: {e}")
            return False

        # Set ceiling override
        data["max_thinking_tokens"] = tokens
        if reason:
            data["ceiling_override_reason"] = reason

        # Save updated metadata
        try:
            os.makedirs(os.path.dirname(meta_path), exist_ok=True)
            with open(meta_path, 'w') as f:
                yaml.dump(data, f, default_flow_style=False, sort_keys=False)

            if verbose:
                msg = f"✅ Ceiling override set: {skill_name} → {tokens} tokens"
                if reason:
                    msg += f". Reason: {reason}"
                print(msg)

            return True

        except Exception as e:
            if verbose:
                print(f"❌ Failed to save ceiling override: {e}")
            return False

    def clear_override(self, skill_name: str, verbose: bool = True) -> bool:
        """
        Clear a ceiling override for a skill.

        Args:
            skill_name: Name of skill
            verbose: Log messages

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

            if "max_thinking_tokens" in data:
                del data["max_thinking_tokens"]
                with open(meta_path, 'w') as f:
                    yaml.dump(data, f, default_flow_style=False, sort_keys=False)

                if verbose:
                    print(f"✅ Ceiling override cleared for {skill_name}")

            if "ceiling_override_reason" in data:
                del data["ceiling_override_reason"]
                with open(meta_path, 'w') as f:
                    yaml.dump(data, f, default_flow_style=False, sort_keys=False)

            return True

        except Exception as e:
            if verbose:
                print(f"❌ Failed to clear ceiling override: {e}")
            return False

    def _load_skill_override(self, skill_name: str) -> Optional[int]:
        """Load max_thinking_tokens override from skill's META.yml."""
        meta_path = self._resolve_meta_path(skill_name)

        if not os.path.exists(meta_path):
            return None

        try:
            with open(meta_path, 'r') as f:
                data = yaml.safe_load(f) or {}

            return data.get("max_thinking_tokens")

        except Exception:
            return None

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

    if len(sys.argv) < 2:
        print("Usage: python ceiling_resolver.py <skill_name> [model] [--set tokens] [--clear]")
        print("Example: python ceiling_resolver.py hyper-architect opus")
        print("Example: python ceiling_resolver.py hyper-architect opus --set 50000")
        sys.exit(1)

    skill_name = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else "opus"

    resolver = CeilingResolver()

    if "--set" in sys.argv:
        idx = sys.argv.index("--set")
        if idx + 1 < len(sys.argv):
            try:
                tokens = int(sys.argv[idx + 1])
                resolver.set_override(skill_name, tokens, verbose=True)
            except ValueError:
                print(f"❌ Invalid token count: {sys.argv[idx + 1]}")
                sys.exit(1)
    elif "--clear" in sys.argv:
        resolver.clear_override(skill_name, verbose=True)
    else:
        ceiling = resolver.get_ceiling(skill_name, model, verbose=True)
        print(f"\nCeiling for {skill_name} ({model}): {ceiling} tokens")
