---
name: new-workflow
description: Converts a desired behavior, prompt, or idea into a properly formatted workflow file in the Hypergraph Coding Agent Framework. Use when the user wants to create a new slash command or workflow.
---

# New Workflow

This skill converts a description of desired agent behavior into a properly formatted workflow file, making it immediately available as a slash command.

## When to use this skill

- When the user asks to create a new "command" or "workflow" from an idea or prompt.
- When the user explicitly runs `/hyper-new-workflow [command_name] [description]`.
- After `/hyper-prompt-engineer` completes and the user chooses to export as a workflow.

## How to use it

1. **Gather Inputs**
   - Extract the command name and its intended behavior from the user's request.
   - If only a description is provided, suggest a short, hyphen-separated name (e.g., `code-review`, `api-design`).
   - Confirm a concise description for the frontmatter.

2. **Draft the Skill File**
   - Create the skill directory: `.agents/skills/<command-name>/`
   - Create the main file: `.agents/skills/<command-name>/SKILL.md`
   - **Mandatory frontmatter** — start the file with:
     ```yaml
     ---
     name: <command-name>
     description: <concise third-person description of what this skill does>
     trigger: /<command-name>
     ---
     ```
   - **Structure the content** using the standard skill format:
     - `## When to use this skill` — trigger conditions
     - `## How to use it` — numbered steps for the agent to follow. Be explicit about agent actions, not just outcomes.

3. **Create IDE Bridge Files**
   - Create a thin bridge in `.claude/commands/<command-name>.md`:
     ```markdown
     ---
     description: "<description matching SKILL.md>"
     ---
     Read `.agents/skills/<command-name>/SKILL.md` and follow its instructions precisely.
     ```
   - Create a thin bridge in `.windsurf/workflows/<command-name>.md`:
     ```markdown
     ---
     description: "<description matching SKILL.md>"
     ---
     Read `.agents/skills/<command-name>/SKILL.md` and follow its instructions precisely.
     ```

4. **Verify and Notify**
   - Confirm all files are well-formed with valid YAML frontmatter.
   - Notify the user: "The `/<command-name>` skill is ready. The SKILL.md is the source of truth in `.agents/skills/<command-name>/`, with IDE bridges in `.claude/commands/` and `.windsurf/workflows/`."
