# Gemini CLI Integration Guide

This document serves as both a manual for human developers and **foundational system instructions** for `gemini-cli`. By reading this file, Gemini CLI understands how to natively orchestrate the Hypergraph Coding Agent Framework.

## 🧑‍💻 For Human Developers

The Hypergraph Framework was originally designed for manual LLM copy-pasting and manual script execution. With **Gemini CLI**, the entire process is automated.

### The Automated Workflow
You do not need to manually open new context windows or run the Python scripts yourself. You simply instruct Gemini CLI to transition through the phases:
1. **Onboarding**: "Gemini, run the Discover phase" or "Run the Baseline phase."
2. **Specification**: "Gemini, activate the Architect skill and let's design a new feature." -> "Now run the Red Team skill on that draft." -> "Let's Resolve the trade-offs."
3. **Execution**: "Gemini, implement MiniPRD-1." Gemini will automatically check the blast radius, write the code, and verify it.

---

## 🤖 For Gemini CLI: System Mandates

**CRITICAL INSTRUCTIONS FOR GEMINI CLI:** When operating in this workspace, you must adhere to the following framework rules:

1. **Skill Integration**: You have native access to the framework's phases via your `activate_skill` tool (`discover`, `baseline`, `architect`, `redteam`, `resolve`, `audit`). When the user requests a phase, activate the corresponding skill immediately and follow its specific procedural guidance strictly.
2. **Autonomous Tool Execution**:
   - You MUST autonomously execute the deterministic python scripts located in `.agents/scripts/` using your `run_shell_command` tool when mandated by a skill's instructions.
   - For example, if a skill requires calculating the blast radius, run `python .agents/scripts/hypergraph_updater.py`.
   - If a skill requires flushing the active workspace, run `python .agents/scripts/archive_specs.py`.
3. **State Management**:
   - Treat `spec/compiled/architecture.yml` as the absolute ground truth for the project's state. 
   - Write temporary drafts and reports to `spec/active/`.
   - Ensure the generated specifications strictly adhere to the templates in `.agents/schemas/`.
4. **Interactive Interviews**: When acting as the Architect, use your `ask_user` tool to gather structured constraints efficiently instead of relying solely on conversational text.
5. **No Probabilistic Traversal**: Do not guess architectural dependencies. Always rely on the outputs of the `hypergraph_updater.py` script to understand the blast radius before executing code modifications.
