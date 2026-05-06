# MiniPRD: CLI Entry Point

**Hypergraph Node ID:** `atom_CLI` (new)
**Parent Node:** `mod_api`

---

## 1. The Confidence Mandate

Before generating any plans or writing code, analyze this document and output a Confidence Score (1-10). If the score is below 9, list strictly the clarifying questions needed to reach 10.

**Prerequisite:** `MiniPRD_MolassemblerAdapter` must be complete and all 6 regression tests must be passing. The CLI wraps the existing `XYZToSMILES` and `OIN3DGenerator` public API — those must be stable before this MiniPRD begins.

---

## 2. Atomic User Stories

* **US-001:** As a pipeline engineer, I want `oin-smiles xyz2oin <path>` to print the OIN string to stdout so that I can chain it with other CLI tools via shell pipes.
* **US-002:** As a pipeline engineer, I want `oin-smiles oin2xyz <oin_string>` to print the XYZ block to stdout so that I can pipe it to a molecular visualizer or further processing.
* **US-003:** As a pipeline engineer, I want non-zero exit codes on any error (file not found, parse error, timeout) so that shell scripts can detect and handle failures.
* **US-004:** As a developer, I want `oin-smiles --help` to describe both subcommands so that users can discover the CLI interface.

---

## 3. Implementation Plan (Task List)

- [ ] Task 1: Create `src/oinsmiles/cli.py` with a `main()` function using `argparse`.
- [ ] Task 2: Set up top-level parser with two subparsers: `xyz2oin` and `oin2xyz`.
  - `xyz2oin`: positional argument `path` (str, path to XYZ file).
  - `oin2xyz`: positional argument `oin` (str, OIN string).
- [ ] Task 3: Implement `xyz2oin` handler:
  1. Check `pathlib.Path(args.path).exists()` → if not, print `FileNotFoundError: {path}` to stderr, `sys.exit(1)`.
  2. Call `XYZToSMILES().convert(args.path)`.
  3. Print result to stdout via `print(oin_string)`.
  4. On any exception: print error to stderr, `sys.exit(1)`.
- [ ] Task 4: Implement `oin2xyz` handler:
  1. Call `OIN3DGenerator().generate(args.oin)`.
  2. Print XYZ block to stdout via `print(xyz_block)`.
  3. On `MolassemblerTimeoutError`: print `Error: Molassembler timed out` to stderr, `sys.exit(2)`.
  4. On any other exception: print error to stderr, `sys.exit(1)`.
- [ ] Task 5: Register CLI entry point in `pyproject.toml` under `[project.scripts]`:
  ```
  oin-smiles = "oinsmiles.cli:main"
  ```
- [ ] Task 6: Run `uv sync` to register the script.
- [ ] Task 7: Verify `oin-smiles --help` works and exits 0.
- [ ] Task 8: Export `main` from `src/oinsmiles/__init__.py` if needed (or leave as internal entry point only — prefer the latter to avoid polluting the public API).

---

## 4. The Negative Space (Constraints)

* **DO NOT** print successful XYZ or OIN output to stderr — stdout only for successful output.
* **DO NOT** add subcommands beyond `xyz2oin` and `oin2xyz` in this MiniPRD.
* **DO NOT** read from stdin by default — file path or OIN string must be explicit positional arguments.
* **DO NOT** add any processing logic to `cli.py` beyond argument parsing and delegation — all logic lives in `XYZToSMILES` and `OIN3DGenerator`.
* **DO NOT** export `main` as part of the public Python API (`from oinsmiles import main`) — it is a CLI entry point only.

---

## 5. Integration Tests & Verification

* **Test 1 (Deterministic):** `oin-smiles --help` → exit code 0; stdout contains "xyz2oin" and "oin2xyz".
* **Test 2 (Deterministic):** `oin-smiles xyz2oin tests/fixtures/cisplatin.xyz` → exit code 0; stdout contains `[Pt_SPL]`.
* **Test 3 (Deterministic):** `oin-smiles xyz2oin /nonexistent_path.xyz` → exit code 1; stderr contains "FileNotFoundError".
* **Test 4 (Deterministic):** `oin-smiles oin2xyz "invalid_oin_string"` → exit code 1; stderr contains error message (does not crash silently).
* **Test 5 (Novel):** `oin-smiles oin2xyz "[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}"` → exit code 0; stdout is a valid XYZ block → [Candidate Artifact routing triggered — XYZ block saved to `tests/candidate_outputs/cli_cisplatin.xyz` for human review]
