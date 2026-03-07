# **Master SOP: Hypergraph Coding Agent Framework Workflow**

**Version:** 1.0 (Strict Sequential Execution Model)

**Core Paradigm:** Spec-First, Deterministic Memory, Automated Auditing

This standard operating procedure dictates the exact sequence of operations for generating, executing, and auditing software utilizing the Hypergraph Coding Agent Framework architecture. Do not deviate from the sequence. Bypassing the specification or auditing phases will result in context collapse, hallucinated requirements, and system graph corruption.

## **Phase \-1: Legacy Onboarding (Brownfield Projects Only)**

If you are implementing this framework in an existing repository, you must map the territory before starting Phase 1\.

1. **Initialize the Hypergraph:** Execute /discover. The agent will scan your code and populate spec/compiled/architecture.yml.  
2. **Review the Map:** Open the YAML file. Verify that the major modules and dependencies were captured correctly.  
3. **Baseline the System:** Execute /baseline. The agent will generate the first SuperPRD.md.  
4. **Verification:** Read the SuperPRD.md. If it accurately reflects what your software does today, you are ready to proceed to **Phase 1: Step 1 (/architect)** to add your first new feature.

## **Phase 0: System Initialization**

### **Greenfield Project Setup**

1. **Clone the Template:** Clone this repository into your target directory. It contains the required .agent/ tooling and empty .gitkeep directories.  
2. **Environment Validation:** Ensure Python 3 is installed. Install the PyYAML dependency required for graph traversal:  
   Bash  
   pip install pyyaml

3. **Execution Permissions:** Ensure the deterministic scripts are executable:  
   Bash  
   chmod \+x .agent/scripts/hypergraph\_updater.py  
   chmod \+x .agent/scripts/archive\_specs.py

## **Phase 1: The Specification Engine**

This phase is entirely conversational. It is designed to extract constraints, map the blast radius of changes, and force architectural trade-offs *before* a single line of code is written.

### **Step 1: Requirements Extraction (/architect)**

1. Open your agentic IDE (Antigravity/Cursor/Claude Code).  
2. Execute the /architect command.  
3. **The Pacing Loop:** The agent will ask a maximum of 2 questions per turn. Answer them thoroughly. Do not attempt to skip phases.  
4. **Artifact:** The agent will generate Draft\_PRD.md in the spec/active/ directory.

### **Step 2: Adversarial Analysis (/redteam)**

1. **CRITICAL:** Open a **New Context Window** (or clear the agent session). Do not let the Red Team agent read the Architect's conversation history.  
2. Execute the /redteam command.  
3. The agent will read Draft\_PRD.md (and architecture.yml if this is an iterative update) to hunt for vulnerabilities, data contract collisions, and missing Non-Functional Requirements (NFRs).  
4. **Artifact:** The agent will generate RedTeam\_Report.md in spec/active/.

### **Step 3: Trade-off Resolution (/resolve)**

1. Open a **New Context Window**.  
2. Execute the /resolve command.  
3. The agent will read the Draft PRD and Red Team Report. It will present you with forced trade-offs (Option A vs. Option B) for every vulnerability found.  
4. **Artifact Generation:** Upon completion, the agent will compile the final SuperPRD.md and the granular MiniPRD\_\[Module\].md files, saving them to spec/compiled/.  
5. **Memory Flush:** The agent will instruct you (or automatically attempt) to execute the archival script. Run it to clear the active workspace:  
   Bash  
   python .agent/scripts/archive\_specs.py \[Feature\_Name\]

## **Phase 2: The Execution Engine (The Build Loop)**

This phase maps to the actual writing of syntax. It is governed by a strict state-machine boundary and deterministic graph traversal.

### **Step 1: The Builder Agent**

1. Open a **New Context Window**.  
2. Prompt your standard Builder agent with a specific MiniPRD target.  
   * *Example Prompt:* "Implement spec/compiled/MiniPRD\_Auth.md. Follow the atomic user stories strictly."  
3. **The Graph Traversal Mandate:** You MUST enforce the following rule on the Builder Agent before it finishes its turn:"Once you have modified the code, execute python .agent/scripts/hypergraph\_updater.py spec/compiled/architecture.yml \[modified\_node\_ids\] to deterministically flag the hypergraph."

### **Step 2: Contract Verification (/audit)**

1. **CRITICAL:** Ensure the Builder Agent has completely stopped generating. *Parallel execution will corrupt the YAML graph.*  
2. Open a **New Context Window**.  
3. Execute the command: /audit spec/compiled/MiniPRD\_\[Target\].md.  
4. The Auditor will evaluate the code against the strict MiniPRD constraints.  
5. **Reconciliation:** If the code passes, the Auditor will semantically rewrite the needs\_review nodes in architecture.yml to match the new inputs/outputs of your codebase, resetting their status to clean.

## **Phase 3: Novel Test Protocol (Human-in-the-Loop)**

Standard CI/CD cannot test subjective or AI-generated outputs. If a MiniPRD contains a Novel Test:

1. The Builder/Auditor will output the result to tests/candidate\_outputs/ and halt execution.  
2. Review the file manually.  
3. If the output is correct, move the file to tests/fixtures/.  
4. Manually update the MiniPRD test definition from novel to deterministic. The system will now automatically run regressions against this approved baseline.

## **Phase 4: Iterative Updates (The Delta Loop)**

When adding new features weeks or months later, the process remains identical.

1. Run /architect. The agent will detect the existing architecture.yml and switch to the **Delta Extraction** protocol.  
2. Run /redteam. It will extract the **Blast Radius** subgraph to identify exactly how your new feature breaks the existing code.  
3. Proceed through /resolve and /audit sequentially.
