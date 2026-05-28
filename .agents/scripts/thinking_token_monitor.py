#!/usr/bin/env python3
"""
Thinking Token Monitor

Resolves thinking token ceiling for skills based on model tier and per-skill overrides.
Logs ceiling application for observability.
"""

import os
import yaml
from typing import Optional


class ThinkingTokenMonitor:
    """Monitors and resolves thinking token ceilings for skills."""

    MODEL_DEFAULTS = {
        "opus": 20000,
        "sonnet": 10000,
        "haiku": 2000,
    }

    def __init__(self, skills_base_dir: str = ".agents/skills"):
        self.skills_base_dir = skills_base_dir

    def monitor(
        self,
        skill_name: str,
        model: str = "opus",
        verbose: bool = True
    ) -> int:
        """
        Monitor and resolve thinking token ceiling for a skill.

        Args:
            skill_name: Name of skill (e.g., "hyper-execute")
            model: Model tier ("opus", "sonnet", "haiku")
            verbose: Log ceiling resolution

        Returns:
            Thinking token ceiling (integer)
        """
        # Try to load skill-specific override from META.yml
        ceiling = self._load_skill_override(skill_name)

        if ceiling is not None:
            if verbose:
                print(f"ℹ️  {skill_name}: Using custom thinking ceiling ({ceiling} tokens)")
            return ceiling

        # Fall back to model tier default
        default = self.MODEL_DEFAULTS.get(model, 2000)

        if verbose:
            print(f"ℹ️  {skill_name}: Using {model} default thinking ceiling ({default} tokens)")

        return default

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
        print("Usage: python thinking_token_monitor.py <skill_name> [model]")
        print("Example: python thinking_token_monitor.py hyper-execute opus")
        sys.exit(1)

    skill_name = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else "opus"

    monitor = ThinkingTokenMonitor()
    ceiling = monitor.monitor(skill_name, model, verbose=True)
    print(f"\nThinking ceiling for {skill_name} ({model}): {ceiling} tokens")
