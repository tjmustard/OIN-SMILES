#!/usr/bin/env python3
"""
Skill Metadata Loader

Loads, validates, and manages skill metadata from META.yml files.
Handles fallbacks for missing metadata and provides schema validation.
"""

import os
import yaml
from dataclasses import dataclass
from typing import Optional, Dict, Any
from datetime import datetime


@dataclass
class SkillMetadata:
    """Represents a skill's metadata configuration."""
    name: str
    assigned_model: str  # "opus", "sonnet", "haiku"
    model_version: Optional[str] = None
    confidence_level: Optional[float] = None
    reasoning_intensity: Optional[str] = None
    max_thinking_tokens: Optional[int] = None
    output_token_budget: Optional[int] = None
    override: Optional[Dict[str, Any]] = None  # {model, set_by, set_at, scope}
    source: str = "META.yml"  # "META.yml", "default", "override"


class MetaLoader:
    """Loads and validates skill metadata from META.yml files."""

    # Default model tier defaults
    TIER_DEFAULTS = {
        "opus": {"max_thinking_tokens": 20000, "output_token_budget": 200000},
        "sonnet": {"max_thinking_tokens": 10000, "output_token_budget": 80000},
        "haiku": {"max_thinking_tokens": 2000, "output_token_budget": 8000},
    }

    # Required fields in META.yml
    REQUIRED_FIELDS = ["name", "assigned_model"]

    def __init__(self, skills_base_dir: str = ".agents/skills"):
        self.skills_base_dir = skills_base_dir

    def load(self, skill_name: str, verbose: bool = True) -> SkillMetadata:
        """
        Load metadata for a skill.

        Args:
            skill_name: Name of skill (e.g., "hyper-execute")
            verbose: Log warnings to console

        Returns:
            SkillMetadata object
        """
        meta_path = self._resolve_meta_path(skill_name)

        if not os.path.exists(meta_path):
            if verbose:
                print(f"⚠️  No META.yml found for {skill_name}. Using default (Opus).")
            return self._default_metadata(skill_name)

        try:
            with open(meta_path, 'r') as f:
                data = yaml.safe_load(f) or {}

            # Validate required fields
            missing = [field for field in self.REQUIRED_FIELDS if field not in data]
            if missing:
                if verbose:
                    print(f"⚠️  META.yml missing required fields: {missing}. Using defaults.")
                return self._default_metadata(skill_name)

            # Build metadata object
            metadata = SkillMetadata(
                name=data.get("name", skill_name),
                assigned_model=data.get("assigned_model", "opus"),
                model_version=data.get("model_version"),
                confidence_level=data.get("confidence_level"),
                reasoning_intensity=data.get("reasoning_intensity"),
                max_thinking_tokens=data.get("max_thinking_tokens"),
                output_token_budget=data.get("output_token_budget"),
                override=data.get("override"),
                source="META.yml"
            )

            # Check for active override
            if metadata.override and metadata.override.get("scope") == "single_run":
                if verbose:
                    print(f"ℹ️  Single-run override active for {skill_name}: "
                          f"{metadata.override.get('model')}")

            return metadata

        except yaml.YAMLError as e:
            if verbose:
                print(f"⚠️  Failed to parse META.yml for {skill_name}: {e}. Using default.")
            return self._default_metadata(skill_name)
        except Exception as e:
            if verbose:
                print(f"⚠️  Error loading metadata for {skill_name}: {e}. Using default.")
            return self._default_metadata(skill_name)

    def _resolve_meta_path(self, skill_name: str) -> str:
        """Resolve path to skill's META.yml file."""
        # Handle both "hyper-execute" and variations
        candidates = [
            os.path.join(self.skills_base_dir, skill_name, "META.yml"),
            os.path.join(self.skills_base_dir, skill_name.replace("-", "_"), "META.yml"),
            os.path.join(self.skills_base_dir, f"hyper-{skill_name}", "META.yml"),
        ]
        return candidates[0]

    def _default_metadata(self, skill_name: str) -> SkillMetadata:
        """Return default metadata (Opus, no version, default ceiling)."""
        return SkillMetadata(
            name=skill_name,
            assigned_model="opus",
            model_version=None,
            confidence_level=None,
            reasoning_intensity=None,
            max_thinking_tokens=self.TIER_DEFAULTS["opus"]["max_thinking_tokens"],
            output_token_budget=self.TIER_DEFAULTS["opus"]["output_token_budget"],
            override=None,
            source="default"
        )

    def validate_schema(self, metadata: Dict[str, Any]) -> tuple:
        """
        Validate metadata against schema.

        Returns:
            (is_valid: bool, errors: List[str])
        """
        errors = []

        # Check required fields
        for field in self.REQUIRED_FIELDS:
            if field not in metadata:
                errors.append(f"Missing required field: {field}")

        # Validate assigned_model
        if "assigned_model" in metadata:
            if metadata["assigned_model"] not in ["opus", "sonnet", "haiku"]:
                errors.append(
                    f"assigned_model must be one of: opus, sonnet, haiku. "
                    f"Got: {metadata['assigned_model']}"
                )

        # Validate optional fields
        if "confidence_level" in metadata:
            if not isinstance(metadata["confidence_level"], (int, float)):
                errors.append("confidence_level must be a number")
            elif not 0 <= metadata["confidence_level"] <= 10:
                errors.append("confidence_level must be between 0 and 10")

        if "reasoning_intensity" in metadata:
            valid_values = ["none", "low", "medium", "high", "very_high"]
            if metadata["reasoning_intensity"] not in valid_values:
                errors.append(
                    f"reasoning_intensity must be one of: {valid_values}. "
                    f"Got: {metadata['reasoning_intensity']}"
                )

        if "max_thinking_tokens" in metadata:
            if not isinstance(metadata["max_thinking_tokens"], int):
                errors.append("max_thinking_tokens must be an integer")
            elif metadata["max_thinking_tokens"] < 0:
                errors.append("max_thinking_tokens must be non-negative")

        return len(errors) == 0, errors


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        skill_name = sys.argv[1]
        loader = MetaLoader()
        metadata = loader.load(skill_name, verbose=True)

        print(f"\nLoaded Metadata for {skill_name}:")
        print(f"  Model: {metadata.assigned_model}")
        print(f"  Version: {metadata.model_version or '(default)'}")
        print(f"  Max Thinking: {metadata.max_thinking_tokens}")
        print(f"  Source: {metadata.source}")
        if metadata.override:
            print(f"  Override: {metadata.override}")
    else:
        print("Usage: python meta_loader.py <skill_name>")
        print("Example: python meta_loader.py hyper-execute")
