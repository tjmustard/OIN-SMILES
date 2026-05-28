#!/usr/bin/env python3
"""
Ceiling Enforcer

Enforces thinking token ceilings by comparing actual consumption against limits.
Detects breaches, logs events, and suggests recovery actions.
"""

from dataclasses import dataclass
from typing import List
from thinking_token_monitor import ThinkingTokenMonitor


@dataclass
class EnforceResult:
    """Result of ceiling enforcement check."""
    skill_name: str
    model: str
    actual_tokens: int
    ceiling: int
    breach: bool
    utilization_percent: float

    def __post_init__(self):
        self.utilization_percent = (self.actual_tokens / self.ceiling * 100) if self.ceiling > 0 else 0


class CeilingEnforcer:
    """Enforces thinking token ceilings and detects breaches."""

    def __init__(self, skills_base_dir: str = ".agents/skills"):
        self.monitor = ThinkingTokenMonitor(skills_base_dir)

    def enforce(
        self,
        skill_name: str,
        model: str,
        actual_thinking_tokens: int,
        verbose: bool = True
    ) -> EnforceResult:
        """
        Enforce thinking token ceiling and detect breaches.

        Args:
            skill_name: Name of skill
            model: Model tier ("opus", "sonnet", "haiku")
            actual_thinking_tokens: Actual tokens consumed
            verbose: Log enforcement result

        Returns:
            EnforceResult with breach status
        """
        ceiling = self.monitor.monitor(skill_name, model, verbose=False)

        result = EnforceResult(
            skill_name=skill_name,
            model=model,
            actual_tokens=actual_thinking_tokens,
            ceiling=ceiling,
            breach=actual_thinking_tokens > ceiling,
            utilization_percent=0
        )
        result.utilization_percent = (actual_thinking_tokens / ceiling * 100) if ceiling > 0 else 0

        if verbose:
            self._log_enforcement(result)

        return result

    def _log_enforcement(self, result: EnforceResult) -> None:
        """Log enforcement result with appropriate formatting."""
        if result.breach:
            self._emit_breach_alert(result)
        else:
            utilization_pct = f"{result.utilization_percent:.0f}%"
            print(f"✅ {result.skill_name}: {result.actual_tokens:,} / {result.ceiling:,} tokens ({utilization_pct})")

    def _emit_breach_alert(self, result: EnforceResult) -> None:
        """Emit detailed breach alert with recovery suggestions."""
        print(f"\n⚠️ Thinking ceiling exceeded: {result.actual_tokens:,} / {result.ceiling:,} tokens ({result.model})")
        print(f"Skill: {result.skill_name}")
        print(f"Utilization: {result.utilization_percent:.1f}%\n")

        suggestions = self._generate_split_suggestions(result)
        print("Suggested split points:")
        for i, suggestion in enumerate(suggestions, 1):
            print(f"  [{i}] {suggestion}")

        print(f"\n💡 Recovery options:")
        print(f"  1. Resubmit Phase 1 with reduced scope")
        print(f"  2. Run: /hyper-config override-ceiling {result.skill_name} {int(result.ceiling * 1.5)} --reason 'Complex task with high reasoning'")
        print(f"  3. Contact framework maintainer if repeatedly hitting ceiling\n")

    def _generate_split_suggestions(self, result: EnforceResult) -> List[str]:
        """Generate suggested split points based on skill."""
        # Generic suggestions; can be customized per skill
        suggestions = [
            "Phase 1: Initial analysis and scoping",
            "Phase 2: Core logic and implementation",
            "Phase 3: Edge cases and validation",
            "Phase 4: Integration and documentation"
        ]
        return suggestions


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 4:
        print("Usage: python ceiling_enforcer.py <skill_name> <model> <actual_tokens>")
        print("Example: python ceiling_enforcer.py hyper-redteam sonnet 12500")
        sys.exit(1)

    skill_name = sys.argv[1]
    model = sys.argv[2]
    try:
        actual_tokens = int(sys.argv[3])
    except ValueError:
        print(f"❌ Invalid token count: {sys.argv[3]}")
        sys.exit(1)

    enforcer = CeilingEnforcer()
    result = enforcer.enforce(skill_name, model, actual_tokens, verbose=True)

    sys.exit(0 if not result.breach else 1)
