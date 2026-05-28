#!/usr/bin/env python3
"""
Token Budget Checker

Enforces output token budgets per MiniPRD to prevent runaway task costs.
Checks estimated vs. ceiling (50k), suggests task splitting, tracks variance.
"""

from dataclasses import dataclass
from typing import List, Optional


@dataclass
class BudgetCheckResult:
    """Result of budget check against 50k ceiling."""
    miniprd_name: str
    estimated_tokens: int
    ceiling: int
    flag: bool
    suggested_subtasks: List[str]
    message: str

    def __post_init__(self):
        if self.estimated_tokens > self.ceiling:
            self.flag = True
            self.message = f"Budget exceeds 50k ceiling: {self.estimated_tokens:,} tokens estimated"
        else:
            self.flag = False
            self.message = f"Budget within limit: {self.estimated_tokens:,} / {self.ceiling:,} tokens"


class BudgetChecker:
    """Checks and enforces output token budgets for MiniPRDs."""

    BUDGET_CEILING = 50000

    def check(
        self,
        miniprd_name: str,
        estimated_output_tokens: int,
        verbose: bool = True
    ) -> BudgetCheckResult:
        """
        Check if estimated output tokens exceed 50k ceiling.

        Args:
            miniprd_name: Name of MiniPRD (e.g., "MiniPRD_TokenBudgetEnforcement")
            estimated_output_tokens: Developer's estimated output token count
            verbose: Log check result

        Returns:
            BudgetCheckResult with flag, message, and split suggestions
        """
        result = BudgetCheckResult(
            miniprd_name=miniprd_name,
            estimated_tokens=estimated_output_tokens,
            ceiling=self.BUDGET_CEILING,
            flag=False,
            suggested_subtasks=[]
        )

        if estimated_output_tokens > self.BUDGET_CEILING:
            result.flag = True
            result.suggested_subtasks = self._generate_split_suggestions(miniprd_name)

        if verbose:
            self._log_check(result)

        return result

    def estimate_tokens(
        self,
        code_lines: int,
        doc_pages: int,
        test_cases: int,
        verbose: bool = True
    ) -> int:
        """
        Estimate output tokens using heuristic formula.

        Formula: code_lines / 40 + doc_pages * 500 + test_cases * 100
        - ~40 tokens per line of code
        - ~500 tokens per documentation page
        - ~100 tokens per test case

        Args:
            code_lines: Estimated lines of code to generate
            doc_pages: Estimated documentation pages
            test_cases: Estimated number of test cases
            verbose: Log estimation

        Returns:
            Estimated output token count
        """
        code_tokens = code_lines / 40
        doc_tokens = doc_pages * 500
        test_tokens = test_cases * 100
        total = int(code_tokens + doc_tokens + test_tokens)

        if verbose:
            print(f"ℹ️  Token estimate: {code_lines} lines @ 40/line + {doc_pages} pages @ 500/page + {test_cases} tests @ 100/test")
            print(f"    = {int(code_tokens)} + {int(doc_tokens)} + {int(test_tokens)} = {total:,} tokens")

        return total

    def _log_check(self, result: BudgetCheckResult) -> None:
        """Log budget check result."""
        status = "❌ FLAG" if result.flag else "✅ PASS"
        print(f"{status}: {result.message}")

        if result.flag and result.suggested_subtasks:
            print(f"\nSuggested splits:")
            for subtask in result.suggested_subtasks:
                print(f"  - {subtask}")

    def _generate_split_suggestions(self, miniprd_name: str) -> List[str]:
        """Generate suggestions for splitting oversized MiniPRDs."""
        return [
            f"[A] Core logic and data structures",
            f"[B] API endpoints and integrations",
            f"[C] Testing and documentation",
        ]


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 4:
        print("Usage: python token_budget_checker.py <code_lines> <doc_pages> <test_cases>")
        print("Example: python token_budget_checker.py 200 5 10")
        print("\nAlternatively:")
        print("Usage: python token_budget_checker.py check <miniprd_name> <estimated_tokens>")
        print("Example: python token_budget_checker.py check MiniPRD_Example 12000")
        sys.exit(1)

    if sys.argv[1] == "check":
        miniprd_name = sys.argv[2]
        try:
            estimated_tokens = int(sys.argv[3])
        except ValueError:
            print(f"❌ Invalid token count: {sys.argv[3]}")
            sys.exit(1)

        checker = BudgetChecker()
        result = checker.check(miniprd_name, estimated_tokens, verbose=True)
        sys.exit(0 if not result.flag else 1)
    else:
        try:
            code_lines = int(sys.argv[1])
            doc_pages = int(sys.argv[2])
            test_cases = int(sys.argv[3])
        except ValueError:
            print("❌ Invalid arguments (must be integers)")
            sys.exit(1)

        checker = BudgetChecker()
        estimated = checker.estimate_tokens(code_lines, doc_pages, test_cases, verbose=True)
        sys.exit(0)
