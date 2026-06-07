---
name: tutorial-generator
description: Generates a markdown tutorial from an integration test or user-provided files through iterative LLM + human collaboration. Tutorials are saved to tutorials/<name>/ with organized supporting files.
trigger: /hyper-tutorial-generator
---

# Tutorial Generator

This skill converts an existing integration test or a set of user-provided files into a polished, step-by-step markdown tutorial. The process is collaborative and interactive — you and the LLM work together to draft and refine each section.

---

## When to use this skill

- After creating an integration test, when you want to document the workflow it demonstrates as a user-facing tutorial.
- When you have example files, screenshots, or code samples that illustrate a framework feature and want to organize them into a structured tutorial.
- When you want to create educational content that walks users through a specific task or workflow.
- Whenever you run `/hyper-tutorial-generator`.

---

## How to use it

### Step 0 — Choose Your Source (HITL Gate #1)

Use **AskUserQuestion** to ask the user which source to use:

**Option A: Existing Integration Test**
- Tutorial will be based on a test file from `tests/integration/`
- You'll select which test and describe what aspects to focus on

**Option B: User-Provided Files**
- Tutorial will be built from files you provide (code, examples, docs, etc.)
- You'll describe the topic and list any supporting files

Wait for the user's selection before proceeding.

---

### Step 1 — Gather Source Material

#### If Integration Test Mode:

1. List all files in `tests/integration/` (numbered).
2. Use **AskUserQuestion** to let the user select which test file to use.
3. Read the selected test file **in full**.
4. Then ask via **AskUserQuestion**: "What aspects of this test should the tutorial focus on? (e.g., the core workflow, specific edge cases, debugging tips)"
5. Wait for the response.

#### If User Files Mode:

1. Use **AskUserQuestion** (two parts):
   - "What is this tutorial about? (1-2 sentences)"
   - "Do you have any input files, code samples, examples, or screenshots to include? If yes, list their file paths or paste their content directly."
2. For any listed files, read them in full (or accept pasted content).
3. Store all provided content in memory for use in later steps.

---

### Step 2 — Propose Tutorial Name and Outline (HITL Gate #2)

Based on the source material, propose:

- **Directory name**: A hyphen-separated name (e.g., `context-clearing-workflow`, `hyper-publish-release-process`, `authentication-integration-test`)
- **Tutorial title**: A human-readable title (e.g., "Context Clearing Workflow", "How to Use Hyper-Publish for Releases")
- **Section outline**: A proposed outline based on the source (e.g., Overview → Prerequisites → Step-by-step → Expected Output → Troubleshooting)

Present the proposal via **AskUserQuestion**:

**Option A: Accept name and outline**
- Proceed with the proposal as-is

**Option B: Edit name / outline**
- Provide your own directory name, title, or reorder/add/remove sections

Wait for the response, then confirm the final name and outline.

---

### Step 3 — Propose Supporting File Structure (HITL Gate #3)

Based on what the source material contains, propose which subdirectories are needed:

- `input_files/` — if the source has example inputs or setup files
- `output_files/` — if there are expected outputs or results to showcase
- `code_samples/` — if there are code snippets or script examples
- `screenshots/` — if visual steps are described (placeholder for user to add images)

Ask via **AskUserQuestion**:

1. "Which supporting subdirectories do you want to include?"
2. "Do you have any additional files to add? (list paths or paste content)"

Based on the response, finalize the subdirectory structure.

---

### Step 4 — Draft Tutorial Section by Section (Iterative HITL Loop)

Work through the outline section by section:

**For each section:**

1. **Draft** the section content based on:
   - The source material (test or user files)
   - The user's stated focus and preferences
   - The template structure in `resources/TUTORIAL_TEMPLATE.md`

2. **Present** the draft to the user using **AskUserQuestion**:
   ```
   Section: <Section Title>
   
   <Drafted content>
   
   Options:
   - Option A: Accept this section
   - Option B: Request changes (describe what you'd like to change)
   - Option C: Regenerate from scratch
   ```

3. **Apply feedback**:
   - If Option A: move to next section
   - If Option B: take the user's feedback and revise the section. Then re-present (go back to step 2 for this section)
   - If Option C: completely regenerate the section with a different approach and re-present

4. **Repeat** until the user accepts the section.

Move on to the next section in the outline.

After all sections are drafted and approved, compile the complete tutorial markdown.

---

### Step 5 — Scaffold Directory and Save Files (HITL Gate #4)

1. Create the full directory structure at project root:
   ```bash
   mkdir -p tutorials/<name>/input_files
   mkdir -p tutorials/<name>/output_files
   mkdir -p tutorials/<name>/code_samples
   # (any other subdirs confirmed in Step 3)
   ```

2. Write the completed tutorial markdown to `tutorials/<name>/tutorial.md`.

3. For any supporting files the user provided:
   - Copy or ask the user to place them in the appropriate subdirectories.
   - If files were pasted as content, write them to the subdirectories now.

---

### Step 6 — Report Result and Next Steps

Display:

```
✅ Tutorial created successfully.

  Location:       tutorials/<name>/
  Main file:      tutorials/<name>/tutorial.md
  Sections:       <N> (e.g., Overview, Prerequisites, Step-by-Step, Troubleshooting)
  Supporting:     input_files/, output_files/, code_samples/
  
  Next steps:
  1. Review the tutorial at tutorials/<name>/tutorial.md
  2. Place any supporting files in the subdirectories
  3. Consider adding this tutorial to docs/index or a README
```

Then ask via **AskUserQuestion**:
- "Would you like to add a summary of this tutorial to `docs/README.md`, create a follow-up tutorial, or finalize?"

---

## Constraints

- **Never overwrite** existing tutorials without explicit user confirmation.
- **Always save** to `tutorials/` at project root (create the directory if it doesn't exist).
- **Always use subdirectories** for supporting files — don't dump everything in the root.
- **Never skip any HITL gate** — the user must approve name, outline, structure, and each section before proceeding.
- **Respect source material** — tutorial should accurately reflect the test or user files it's based on; don't add unsupported claims.

---

## Template Reference

Use `.agents/skills/hyper-tutorial-generator/resources/TUTORIAL_TEMPLATE.md` as a reference when drafting each section. The standard sections are:
- **Overview** — what this tutorial demonstrates and why it matters
- **Prerequisites** — what the user needs before starting
- **Step-by-Step** — numbered, detailed steps with explanations
- **Expected Output** — what the user should see when done
- **Troubleshooting** — common problems and solutions
- **Related** — links to related skills, docs, or tutorials
