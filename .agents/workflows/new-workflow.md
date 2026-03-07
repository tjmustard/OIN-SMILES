---
description: Create new slash commands on the fly.
---

# /new-workflow Workflow

Follow these steps to predictably convert a desired behavior, prompt, or existing workflow pattern into a properly formatted workflow file for the Hypergraph Coding Agent Framework.

## When to use this workflow
- When the user asks you to create a "command" or "workflow" out of an idea or prompt.
- When the user explicitly runs the `/newcommand [command_name] [description]` slash command.

## How to use it

1. **Gather Inputs**: 
   - Extract the command name (`[command_name]`) and its description or intended behavior.
   - If only a description is provided, suggest a short, hyphen-separated name for the workflow (e.g., `code-review` or `api-design`).
   - Determine a concise description of what the workflow does to include in its frontmatter.

2. **Draft the Workflow File**:
   - Create a new markdown file named after the command in the `.agents/workflows/` directory: `.agents/workflows/<command-name>.md`.
   - **Mandatory Frontmatter**: Start the file with the following YAML frontmatter:
     ```yaml
     ---
     description: <short title, e.g. how to deploy the application>
     ---
     ```
   - **Structure the Content**: Write out the specific steps on how to execute this workflow using a numbered list or bullet points. Include clear, actionable agentic instructions describing exactly what the agent should do when the `/command-name` is invoked.

3. **Verify and Notify**:
   - Confirm the file is well-formed, valid YAML is used in the frontmatter, and it has been written to the `.agents/workflows/` directory.
   - Notify the user that the new command is ready to use and summarize its intended behavior.