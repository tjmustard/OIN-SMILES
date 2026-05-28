#!/usr/bin/env python3
"""
Token Budget Reconciler

Reconciles estimated vs. actual output tokens post-execution.
Tracks variance and flags estimation model drift (>20% variance).
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class ReconciliationResult:
    """Result of post-execution token budget reconciliation."""
    miniprd_name: str
    estimated_tokens: int
    actual_tokens: int
    variance_percent: float
    investigation_required: bool
    message: str

    def __post_init__(self):
        if self.estimated_tokens > 0:
            self.variance_percent = ((self.actual_tokens - self.estimated_tokens) / self.estimated_tokens * 100)
        else:
            self.variance_percent = 0

        self.investigation_required = abs(self.variance_percent) > 20


class TokenBudgetReconciler:
    """Reconciles estimated vs. actual output tokens."""

    VARIANCE_THRESHOLD = 20  # percent

    def reconcile(
        self,
        miniprd_name: str,
        estimated_tokens: int,
        actual_tokens: int,
        verbose: bool = True
    ) -> ReconciliationResult:
        """
        Reconcile estimated vs. actual output tokens.

        Args:
            miniprd_name: Name of MiniPRD
            estimated_tokens: Developer's pre-execution estimate
            actual_tokens: Actual tokens consumed
            verbose: Log reconciliation result

        Returns:
            ReconciliationResult with variance metrics
        """
        result = ReconciliationResult(
            miniprd_name=miniprd_name,
            estimated_tokens=estimated_tokens,
            actual_tokens=actual_tokens,
            variance_percent=0,
            investigation_required=False
        )

        if verbose:
            self._log_reconciliation(result)

        return result

    def _log_reconciliation(self, result: ReconciliationResult) -> None:
        """Log reconciliation result with variance analysis."""
        direction = "+" if result.variance_percent > 0 else ""
        variance_str = f"{direction}{result.variance_percent:.1f}%"

        print(f"📊 {result.miniprd_name}: {result.actual_tokens:,} / {result.estimated_tokens:,} tokens (variance: {variance_str})")

        if result.investigation_required:
            print(f"⚠️  Variance exceeds {self.VARIANCE_THRESHOLD}% threshold — investigate estimation model accuracy")
            print(f"   • Review estimation formula (code_lines/40 + doc_pages*500 + test_cases*100)")
            print(f"   • Update heuristic if pattern emerges across multiple MiniPRDs")
        else:
            print(f"✅ Variance within acceptable range ({self.VARIANCE_THRESHOLD}%)")


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 4:
        print("Usage: python token_budget_reconciler.py <miniprd_name> <estimated> <actual>")
        print("Example: python token_budget_reconciler.py MiniPRD_Example 8000 7200")
        sys.exit(1)

    miniprd_name = sys.argv[1]
    try:
        estimated = int(sys.argv[2])
        actual = int(sys.argv[3])
    except ValueError:
        print("❌ Invalid token counts (must be integers)")
        sys.exit(1)

    reconciler = TokenBudgetReconciler()
    result = reconciler.reconcile(miniprd_name, estimated, actual, verbose=True)

    sys.exit(0 if not result.investigation_required else 1)
