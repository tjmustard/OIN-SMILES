#!/usr/bin/env python3
"""
Model Routing Matrix

Routes skills to appropriate Claude models (Opus/Sonnet/Haiku) based on META.yml
metadata and heuristic rules. Handles fallback, version pinning, and user overrides.
"""

import os
import yaml
from dataclasses import dataclass
from typing import Optional, Dict, Any, List
from enum import Enum


class ModelTier(Enum):
    """Available Claude model tiers."""
    HAIKU = "haiku"
    SONNET = "sonnet"
    OPUS = "opus"


@dataclass
class ModelRouterResponse:
    """Response from model routing decision."""
    assigned_model: str  # "opus", "sonnet", "haiku"
    model_version: str  # Full model ID: "claude-opus-4-7"
    max_thinking_tokens: int
    source: str  # "META.yml", "heuristic_suggestion", "fallback", "user_override"
    confidence: float = 1.0
    warning_message: Optional[str] = None


class ModelRouter:
    """Routes skills to appropriate Claude models."""

    def __init__(self, rules_path: str = ".agents/config/heuristic_rules.yaml"):
        self.rules_path = rules_path
        self.rules = self._load_rules()
        self.model_mapping = self.rules.get("model_mapping", {})
        self.fallback_policy = self.rules.get("fallback", {})
        self.versions = self.rules.get("versions", {})

    def _load_rules(self) -> Dict[str, Any]:
        """Load heuristic rules from YAML configuration."""
        if not os.path.exists(self.rules_path):
            return {"model_mapping": {}, "fallback": {}, "versions": {}}

        try:
            with open(self.rules_path, 'r') as f:
                return yaml.safe_load(f) or {}
        except Exception as e:
            print(f"Warning: Failed to load heuristic rules: {e}")
            return {"model_mapping": {}, "fallback": {}, "versions": {}}

    def route(
        self,
        skill_name: str,
        user_override: Optional[str] = None,
        verbose: bool = True
    ) -> ModelRouterResponse:
        """
        Route a skill to an appropriate Claude model.

        Args:
            skill_name: Name of skill (e.g., "hyper-execute", "hyper-architect")
            user_override: Temporary user-specified model tier ("opus", "sonnet", "haiku")
            verbose: Log warnings and alerts

        Returns:
            ModelRouterResponse with assigned model and metadata
        """
        # Handle user override
        if user_override:
            return self._handle_user_override(user_override)

        # Try to read META.yml
        meta_path = self._resolve_meta_path(skill_name)

        if os.path.exists(meta_path):
            return self._route_from_metadata(meta_path, verbose)
        else:
            return self._handle_missing_metadata(skill_name, verbose)

    def _resolve_meta_path(self, skill_name: str) -> str:
        """Resolve the path to a skill's META.yml file."""
        # Handle both "hyper-execute" and "hyper-execute" formats
        skill_dir = skill_name.replace("-", "_")
        candidates = [
            f".agents/skills/{skill_name}/META.yml",
            f".agents/skills/{skill_dir}/META.yml",
            f".agents/skills/hyper-{skill_name}/META.yml",
        ]
        return candidates[0]  # Return first candidate; caller checks existence

    def _route_from_metadata(self, meta_path: str, verbose: bool) -> ModelRouterResponse:
        """Route based on skill's META.yml metadata."""
        try:
            with open(meta_path, 'r') as f:
                metadata = yaml.safe_load(f) or {}

            assigned_model = metadata.get("assigned_model", "opus").lower()
            model_version = metadata.get(
                "model_version",
                self.model_mapping.get(assigned_model, {}).get("model_id", "")
            )
            max_thinking = metadata.get(
                "max_thinking_tokens",
                self.model_mapping.get(assigned_model, {}).get("max_thinking_tokens", 1024)
            )

            # Check version compatibility
            warning = self._check_version_compatibility(assigned_model, model_version)
            if warning and verbose:
                print(warning)

            return ModelRouterResponse(
                assigned_model=assigned_model,
                model_version=model_version,
                max_thinking_tokens=max_thinking,
                source="META.yml",
                warning_message=warning
            )
        except Exception as e:
            if verbose:
                print(f"Warning: Failed to read {meta_path}: {e}")
            return self._handle_missing_metadata(os.path.dirname(meta_path), verbose)

    def _handle_missing_metadata(self, skill_name: str, verbose: bool) -> ModelRouterResponse:
        """Handle fallback when META.yml is missing."""
        default_model = self.fallback_policy.get("default_model", "opus")
        message = self.fallback_policy.get(
            "alert_message",
            f"⚠️ {skill_name} has no model assignment. Routing to Opus (fallback)."
        ).format(skill=skill_name)

        if verbose and self.fallback_policy.get("alert_terminal", True):
            print(message)

        model_version = self.model_mapping.get(default_model, {}).get("model_id", "")
        max_thinking = self.model_mapping.get(default_model, {}).get("max_thinking_tokens", 1024)

        return ModelRouterResponse(
            assigned_model=default_model,
            model_version=model_version,
            max_thinking_tokens=max_thinking,
            source="fallback",
            warning_message=message
        )

    def _handle_user_override(self, override_tier: str) -> ModelRouterResponse:
        """Handle user-specified temporary model override."""
        tier = override_tier.lower()
        if tier not in self.model_mapping:
            return self._handle_missing_metadata(f"invalid_override[{override_tier}]", True)

        model_version = self.model_mapping[tier].get("model_id", "")
        max_thinking = self.model_mapping[tier].get("max_thinking_tokens", 1024)

        return ModelRouterResponse(
            assigned_model=tier,
            model_version=model_version,
            max_thinking_tokens=max_thinking,
            source="user_override",
            warning_message="User override applied (single run). Reverting to META.yml assignment after execution."
        )

    def _check_version_compatibility(
        self,
        model_tier: str,
        model_version: str
    ) -> Optional[str]:
        """Check if model version is supported; warn if retired."""
        tier_versions = self.versions.get(model_tier, {})
        supported = tier_versions.get("supported", [])
        retired = tier_versions.get("retired", [])

        if model_version in retired:
            return (
                f"⚠️ Model version '{model_version}' is retired. "
                f"Consider updating META.yml to '{tier_versions.get('current', '')}'"
            )

        if supported and model_version not in supported:
            return (
                f"⚠️ Model version '{model_version}' not in supported list. "
                f"May be unavailable or have different behavior."
            )

        return None


class SkillMetadataGenerator:
    """Generates META.yml template for skills."""

    @staticmethod
    def create_template(
        skill_name: str,
        assigned_model: str = "opus",
        description: str = ""
    ) -> str:
        """Generate META.yml template for a skill."""
        return f"""---
name: {skill_name}
description: {description}

# Model assignment: haiku | sonnet | opus
assigned_model: {assigned_model}

# Optional: pin to specific model version
# If omitted, uses current version from heuristic_rules.yaml
# model_version: "claude-{assigned_model}-4-5-20251001"

# Optional: override max thinking tokens
# max_thinking_tokens: 8192

# Optional: metadata for heuristic classifier
# tool_calls: 5  # Expected number of tool calls
# reasoning_intensity: "medium"  # none | low | medium | high | very_high
# output_determinism: "high"  # high | medium | low
"""


if __name__ == "__main__":
    import sys

    router = ModelRouter()

    if len(sys.argv) > 1:
        skill_name = sys.argv[1]
        override = sys.argv[2] if len(sys.argv) > 2 else None

        result = router.route(skill_name, override, verbose=True)

        print(f"\nRouting Result:")
        print(f"  Skill: {skill_name}")
        print(f"  Model: {result.assigned_model}")
        print(f"  Version: {result.model_version}")
        print(f"  Max Thinking: {result.max_thinking_tokens}")
        print(f"  Source: {result.source}")
        print(f"  Confidence: {result.confidence}")
    else:
        print("Usage: python model_router.py <skill_name> [user_override_tier]")
        print("Example: python model_router.py hyper-architect")
        print("Example: python model_router.py hyper-execute sonnet")
