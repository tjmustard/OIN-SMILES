---
name: new-workflow
description: Converts a desired behavior, prompt, or idea into a properly formatted workflow file in the Hypergraph Coding Agent Framework. Use when the user wants to create a new slash command or workflow.
---

# New Workflow

This skill converts a description of desired agent behavior into a properly formatted workflow file, making it immediately available as a slash command.

## When to use this skill

- When the user asks to create a new "command" or "workflow" from an idea or prompt.
- When the user explicitly runs `/new-workflow [command_name] [description]`.
- After `/prompt-engineer` completes and the user chooses to export as a workflow.

## How to use it

1. **Gather Inputs**
   - Extract the command name and its intended behavior from the user's request.
   - If only a description is provided, suggest a short, hyphen-separated name (e.g., `code-review`, `api-design`).
   - Confirm a concise description for the frontmatter.

2. **Draft the Workflow File**
   - Create a new markdown file at `.agents/workflows/<command-name>.md`.
   - **Mandatory frontmatter** — start the file with:
     ```yaml
     ---
     description: <short description of what this workflow does>
     ---
     ```
   - **Structure the content**: Write specific, actionable steps for what the agent should do when the command is invoked. Use numbered lists or bullet points. Be explicit about agent actions, not just outcomes.

3. **Verify and Notify**
   - Confirm the file is well-formed with valid YAML frontmatter.
   - Confirm it has been written to `.agents/workflows/`.
   - Notify the user: "The `/<command-name>` workflow is ready. Here's what it does: [summary]."
