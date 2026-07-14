#!/usr/bin/env python3
"""
Calculate lines of code and characters of code per folder, separated by file type.
Provides a summary of different categories of code.
"""

import os
from collections import defaultdict
from pathlib import Path

# Configuration
ROOT_DIR = Path(__file__).resolve().parent.parent

# Which extensions to track.
TRACKED_EXTENSIONS = {
    ".py",
    ".md",
    ".sh",
    ".yml",
    ".yaml",
    ".toml",
    ".json",
    ".cff",
    ".txt",
    ".csv",
}

# Directories to exclude from the traversal entirely.
EXCLUDED_DIRS = {
    ".git",
    ".venv",
    ".pytest_cache",
    ".ruff_cache",
    "build",
    ".claude",
    ".github",
    "__pycache__",
    "tmCAT-tmPHOTO_xyz_dataset",
}


def get_metrics():
    # folder -> extension -> {'lines': 0, 'chars': 0}
    metrics = defaultdict(lambda: defaultdict(lambda: {"lines": 0, "chars": 0}))

    for root, dirs, files in os.walk(ROOT_DIR):
        # Exclude specified directories and hidden directories
        dirs[:] = [d for d in dirs if d not in EXCLUDED_DIRS and not d.startswith(".")]

        rel_dir = os.path.relpath(root, ROOT_DIR)
        if rel_dir == ".":
            rel_dir = "root"

        for file in files:
            file_path = Path(root) / file
            ext = file_path.suffix.lower()

            # Use filename as extension for extensionless files we care about (like LICENSE)
            if not ext and file_path.name in {"LICENSE", "Makefile"}:
                ext = file_path.name

            if ext in TRACKED_EXTENSIONS or ext in {"LICENSE", "Makefile"}:
                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        content = f.read()
                        chars = len(content)
                        lines = len(content.splitlines())

                        if chars > 0 or lines > 0:
                            metrics[rel_dir][ext]["lines"] += lines
                            metrics[rel_dir][ext]["chars"] += chars
                except Exception:
                    # Skip files that can't be read as utf-8 text
                    pass

    return metrics


def main():
    metrics = get_metrics()

    print(f"{'Folder':<35} | {'Type':<10} | {'Lines':<10} | {'Characters'}")
    print("-" * 75)

    total_lines = defaultdict(int)
    total_chars = defaultdict(int)

    for folder in sorted(metrics.keys()):
        for ext in sorted(metrics[folder].keys()):
            lines = metrics[folder][ext]["lines"]
            chars = metrics[folder][ext]["chars"]
            total_lines[ext] += lines
            total_chars[ext] += chars
            print(f"{folder:<35} | {ext:<10} | {lines:<10} | {chars}")

    print("-" * 75)
    print("TOTALS BY FILE TYPE:")
    for ext in sorted(total_lines.keys()):
        print(f"  {ext:<33} | TYPE TOTAL | {total_lines[ext]:<10} | {total_chars[ext]}")

    print("\n" + "=" * 75)
    print("CATEGORY SUMMARY (Lines of Code)")
    print("=" * 75)

    categories = {
        "Actual Code (src/)": lambda f: f.startswith("src"),
        "Testing Code (tests/)": lambda f: f.startswith("tests"),
        "Tools Code (tools/)": lambda f: f.startswith("tools"),
        "Examples Code (examples/)": lambda f: f.startswith("examples"),
        "Devtools (devtools/)": lambda f: f.startswith("devtools"),
        "Documentation (docs/)": lambda f: f.startswith("docs"),
        "Root & Other": lambda f: (
            not any(
                f.startswith(p) for p in ["src", "tests", "tools", "examples", "devtools", "docs"]
            )
        ),
    }

    for cat_name, condition in categories.items():
        cat_metrics = defaultdict(int)
        for folder, ext_data in metrics.items():
            if condition(folder):
                for ext, data in ext_data.items():
                    cat_metrics[ext] += data["lines"]

        if cat_metrics:
            summary_str = ", ".join(
                f"{v} {k}"
                for k, v in sorted(cat_metrics.items(), key=lambda item: item[1], reverse=True)
            )
            print(f"{cat_name:<25} : {summary_str}")


if __name__ == "__main__":
    main()
