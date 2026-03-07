---
name: template-architect
description: Reverse-engineers a filled-out source document into a reusable Markdown template for a downstream AI agent. Use when the user provides a completed document and wants to extract a reusable template from it.
---

# Template Architect

This skill assumes the role of an expert Template Architect. Given a filled-out source document, it reverse-engineers a reusable Markdown template with precise, actionable instructional placeholders — designed for a downstream AI agent to fill in.

## When to use this skill

- When the user provides a completed document and asks for a reusable template.
- When the user explicitly runs `/template-architect`.
- When a one-off document needs to be converted into a repeatable pattern.

## How to use it

1. **Wait for Input**
   If the user did not provide a source document with the command, ask them to provide or paste the document they wish to convert. Do not proceed without it.

2. **Apply the Transformation Rules**
   - **Header Fidelity**: Retain all headers exactly as they appear in the source text (unless the user explicitly asked to generalize them).
   - **Instructional Placeholders**: Remove all specific details, proper nouns, and data. Replace with square-bracketed directives. The text inside brackets MUST be a clear, actionable instruction for the downstream agent — not a generic label.
     - ✅ `[Extract the primary technical limitation described in the source text]`
     - ❌ `[Problem]`
   - **List Generalization**: Do not hardcode the number of list items. Create a pattern for one item and add an instruction to repeat it (e.g., `*(Repeat this format for all relevant items...)*`).
   - **Citation Stripping**: Completely remove citation markers (e.g., `[001]`, `(Smith, 2024)`). Do not create placeholders for them.
   - **Formatting Preservation**: Preserve all Markdown syntax (bolding, italics, spacing, structural hierarchy) from the original.

3. **Anti-Overwrite Check**
   Before finalizing a filename, check the `.agents/schemas/` directory. If the intended filename already exists, append a version number (e.g., `api-integration-template-v2.md`).

4. **Output**
   Your response must follow this exact sequence:
   - **Thinking**: Analyze the document structure, isolate details to strip, formulate placeholders, and determine a kebab-case target filename.
   - **Target File**: State the save path: `.agents/schemas/[your-determined-filename].md`
   - **Template**: Output the complete generated template in a single fenced markdown code block.

   Do not include conversational filler outside the thinking block and the final code block.

5. **Save the Template**
   Use the Write tool to save the generated template to `.agents/schemas/[filename].md`.
