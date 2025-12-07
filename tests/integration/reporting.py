
import sys
from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class TestResult:
    name: str
    status: str # "PASS" or "FAIL"
    details: str = ""
    expected: Optional[str] = None
    got: Optional[str] = None

class VerificationReporter:
    def __init__(self, title: str = "Verification Report"):
        self.title = title
        self.results: List[TestResult] = []

    def log_success(self, name: str, details: str = ""):
        self.results.append(TestResult(name, "PASS", details))

    def log_failure(self, name: str, reason: str, expected: Optional[str] = None, got: Optional[str] = None):
        self.results.append(TestResult(name, "FAIL", reason, expected, got))

    def print_summary(self):
        print(f"\n# {self.title}\n")
        print("| Test Case | Status | Details |")
        print("| :--- | :--- | :--- |")
        
        passed = 0
        failed = 0
        
        for res in self.results:
            status_icon = "✅ PASS" if res.status == "PASS" else "❌ FAIL"
            if res.status == "PASS":
                passed += 1
            else:
                failed += 1
            
            # Format Details
            details = res.details
            if res.expected or res.got:
                details += "<br>"
                if res.expected: details += f"**Expected:** `{res.expected}`<br>"
                if res.got: details += f"**Got:** `{res.got}`"
            
            print(f"| {res.name} | {status_icon} | {details} |")
            
        print(f"\n**Total:** {passed} Passed, {failed} Failed.")

        if failed > 0:
            print("\n[FAILURE] Test suite failed.")
            # We don't exit(1) here to allow calling script to handle exit code if needed, 
            # but usually this is the end.
