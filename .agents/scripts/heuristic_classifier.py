#!/usr/bin/env python3
"""
Heuristic Classifier

Analyzes skill characteristics and suggests appropriate Claude model tier.
Helps developers assign correct model without guessing.
"""

import yaml
import re
from dataclasses import dataclass
from typing import Optional, List, Dict, Any
from enum import Enum


class ModelTier(Enum):
    """Available model tiers."""
    HAIKU = "haiku"
    SONNET = "sonnet"
    OPUS = "opus"


@dataclass
class ClassifierSuggestion:
    """Suggestion for appropriate model tier."""
    suggested_tier: str  # "haiku", "sonnet", "opus"
    confidence: float  # 0.0 to 1.0
    reasoning: List[str]  # Human-readable reasons for suggestion
    alternatives: Dict[str, float]  # Confidence for other tiers
    rule_matched: str  # Which rule matched best


class SkillSpec:
    """Represents a skill specification for analysis."""

    def __init__(
        self,
        name: str,
        description: str,
        tool_calls: Optional[int] = None,
        reasoning_intensity: Optional[str] = None,
        output_determinism: Optional[str] = None,
        content: Optional[str] = None
    ):
        self.name = name
        self.description = description
        self.tool_calls = tool_calls
        self.reasoning_intensity = reasoning_intensity
        self.output_determinism = output_determinism
        self.content = content or ""

    def infer_characteristics(self) -> None:
        """Infer characteristics from skill content if not provided."""
        if not self.content:
            return

        # Count tool calls mentioned in content
        tool_pattern = r"(Read|Edit|Write|Bash|Agent)\("
        tool_count = len(re.findall(tool_pattern, self.content))
        if not self.tool_calls:
            self.tool_calls = tool_count

        # Infer reasoning intensity from keywords
        heavy_keywords = ["complex", "adversarial", "architecture", "trade", "red team", "strategy"]
        medium_keywords = ["analyze", "generate", "optimize", "coordinate", "decide"]
        light_keywords = ["format", "extract", "simple", "read", "list"]

        content_lower = (self.description + " " + self.content).lower()

        if not self.reasoning_intensity:
            if any(kw in content_lower for kw in heavy_keywords):
                self.reasoning_intensity = "high"
            elif any(kw in content_lower for kw in medium_keywords):
                self.reasoning_intensity = "medium"
            else:
                self.reasoning_intensity = "low"

        # Infer output determinism
        if not self.output_determinism:
            deterministic_keywords = ["format", "extract", "parse", "list", "read"]
            if any(kw in content_lower for kw in deterministic_keywords):
                self.output_determinism = "high"
            else:
                self.output_determinism = "medium"


class HeuristicClassifier:
    """Classifies skills into model tiers using heuristic rules."""

    def __init__(self, rules_path: str = ".agents/config/heuristic_rules.yaml"):
        self.rules_path = rules_path
        self.rules = self._load_rules()
        self.tiers = self.rules.get("tiers", {})

    def _load_rules(self) -> Dict[str, Any]:
        """Load heuristic rules from YAML."""
        try:
            with open(self.rules_path, 'r') as f:
                return yaml.safe_load(f) or {}
        except FileNotFoundError:
            print(f"Warning: Rules file not found: {self.rules_path}")
            return {"tiers": {}}

    def suggest(self, skill_spec: SkillSpec) -> ClassifierSuggestion:
        """
        Suggest appropriate model tier for a skill.

        Args:
            skill_spec: SkillSpec with skill characteristics

        Returns:
            ClassifierSuggestion with tier recommendation and confidence
        """
        # Infer missing characteristics from content
        skill_spec.infer_characteristics()

        # Score each tier
        scores = {}
        matched_rules = {}

        for tier_name in ["haiku", "sonnet", "opus"]:
            score, rule = self._score_tier(skill_spec, tier_name)
            scores[tier_name] = score
            matched_rules[tier_name] = rule

        # Find best match
        best_tier = max(scores, key=scores.get)
        best_score = scores[best_tier]

        # Build reasoning
        reasoning = [
            f"Tool calls: {skill_spec.tool_calls or 'unknown'}",
            f"Reasoning intensity: {skill_spec.reasoning_intensity or 'unknown'}",
            f"Output determinism: {skill_spec.output_determinism or 'unknown'}",
            f"Matched rule: {matched_rules[best_tier]}",
        ]

        return ClassifierSuggestion(
            suggested_tier=best_tier,
            confidence=min(best_score, 1.0),
            reasoning=reasoning,
            alternatives={
                tier: min(score, 1.0)
                for tier, score in scores.items()
                if tier != best_tier
            },
            rule_matched=matched_rules[best_tier]
        )

    def _score_tier(self, spec: SkillSpec, tier_name: str) -> tuple:
        """Score how well a skill matches a tier. Returns (score, matched_rule)."""
        tier_rules = self.tiers.get(tier_name, {}).get("conditions", [])

        best_match_score = 0.0
        best_match_rule = "no_rule_matched"

        for rule in tier_rules:
            score = self._match_rule(spec, rule)
            if score > best_match_score:
                best_match_score = score
                best_match_rule = rule.get("name", "unnamed_rule")

        return best_match_score, best_match_rule

    def _match_rule(self, spec: SkillSpec, rule: Dict[str, Any]) -> float:
        """Score how well a skill matches a rule (0.0 to 1.0)."""
        score = 0.0
        matched_conditions = 0
        total_conditions = 0

        # Check tool_calls range
        if "tool_calls" in rule:
            total_conditions += 1
            min_calls, max_calls = rule["tool_calls"]
            if spec.tool_calls is not None:
                if min_calls <= spec.tool_calls <= max_calls:
                    score += rule.get("confidence", 0.9)
                    matched_conditions += 1

        # Check reasoning_intensity
        if "reasoning_intensity" in rule:
            total_conditions += 1
            allowed = rule["reasoning_intensity"]
            if not isinstance(allowed, list):
                allowed = [allowed]
            if spec.reasoning_intensity in allowed:
                score += rule.get("confidence", 0.9)
                matched_conditions += 1

        # Check output_determinism
        if "output_determinism" in rule:
            total_conditions += 1
            allowed = rule["output_determinism"]
            if not isinstance(allowed, list):
                allowed = [allowed]
            if spec.output_determinism in allowed:
                score += rule.get("confidence", 0.9)
                matched_conditions += 1

        # Check tools_allowed
        if "tools_allowed" in rule:
            total_conditions += 1
            # Simplified: just check if any allowed tools are mentioned
            content = (spec.description + " " + spec.content).lower()
            allowed_tools = [t.lower() for t in rule["tools_allowed"]]
            if any(tool in content for tool in allowed_tools):
                matched_conditions += 1

        # Check keywords
        if "keywords" in rule:
            total_conditions += 1
            content = (spec.description + " " + spec.content).lower()
            keywords = [kw.lower() for kw in rule["keywords"]]
            if any(kw in content for kw in keywords):
                score += rule.get("confidence", 0.9)
                matched_conditions += 1

        # Normalize score
        if total_conditions > 0:
            return score / total_conditions
        else:
            return rule.get("confidence", 0.5)


if __name__ == "__main__":
    import sys

    classifier = HeuristicClassifier()

    if len(sys.argv) > 1:
        skill_name = sys.argv[1]
        description = sys.argv[2] if len(sys.argv) > 2 else ""

        spec = SkillSpec(
            name=skill_name,
            description=description
        )
        spec.infer_characteristics()

        suggestion = classifier.suggest(spec)

        print(f"\nClassification Result:")
        print(f"  Skill: {skill_name}")
        print(f"  Suggested Tier: {suggestion.suggested_tier}")
        print(f"  Confidence: {suggestion.confidence:.2f}")
        print(f"  Matched Rule: {suggestion.rule_matched}")
        print(f"\nReasoning:")
        for reason in suggestion.reasoning:
            print(f"  - {reason}")
        print(f"\nAlternatives:")
        for tier, conf in suggestion.alternatives.items():
            print(f"  - {tier}: {conf:.2f}")
    else:
        print("Usage: python heuristic_classifier.py <skill_name> [description]")
        print("Example: python heuristic_classifier.py hyper-architect 'Requirements interview'")
