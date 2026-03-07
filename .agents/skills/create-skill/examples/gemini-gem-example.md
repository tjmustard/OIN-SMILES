# Example Input:
Create a skill called "commit-formatter" that takes a messy git diff and turns it into a structured semantic commit message. Use the conventional commits format. It should be used whenever I ask you to write a commit message.

Make sure to include a sample regex for extracting scopes in `scripts/extract_scope.py`.

# Example Output `SKILL.md` (in `.agents/skills/commit-formatter/SKILL.md`):
---
name: commit-formatter
description: Generates a structured semantic commit message based on a provided git diff using Conventional Commits. Use whenever asked to write or format a commit message.
---

# Commit Formatter

This skill helps standardizes commit messages across the repository using the Conventional Commits specification.

## When to use this skill
- When the user provides a `git diff` and asks for a commit message.
- When the user explicitly asks to format a commit message.

## How to use it
1. **Analyze the Diff**: Review the provided diff to understand the core changes.
2. **Determine the Type**: Categorize the change into one of the following: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `chore`.
3. **Determine the Scope**: Identify the module or component affected. You can use `scripts/extract_scope.py` for hints.
4. **Draft the Subject**: Write a concise subject line in the imperative mood, all lowercase, without a period at the end.
5. **Draft the Body**: Provide a detailed description of *why* the change was made, detailing the problem and the solution.
6. **Format Output**: Present the final message in a code block for the user to easily copy.


# Example Output Auxiliary File (in `.agents/skills/commit-formatter/scripts/extract_scope.py`)
```python
import re
import sys

def extract_scope(diff_text):
    # Example logic to find folder names in diff
    matches = re.findall(r'a/([^/]+)/', diff_text)
    if matches:
        return matches[0]
    return "core"

if __name__ == "__main__":
    diff = sys.stdin.read()
    print(extract_scope(diff))
```
