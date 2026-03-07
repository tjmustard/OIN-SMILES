# MiniPRD: Molassembler Spike

**Hypergraph Node ID:** `atom_MolassemblerSpike` (temporary spike — not added to architecture.yml)
**Parent Node:** `mod_generation`

---

## 1. The Confidence Mandate

Before generating any plans or writing code, analyze this document and output a Confidence Score (1-10). If the score is below 9, list strictly the clarifying questions needed to reach 10.

**Key uncertainty to resolve:** Whether `scine.molassembler` Python bindings are picklable for use with `concurrent.futures.ProcessPoolExecutor`. This MiniPRD must complete before `MiniPRD_MolassemblerAdapter` begins.

---

## 2. Atomic User Stories

* **US-001:** As a developer, I want to confirm `scine.molassembler` imports correctly in the project environment so that I know the SCINE stack is properly installed.
* **US-002:** As a developer, I want to confirm that the module-level worker function for `ProcessPoolExecutor` is picklable so that the GIL-safe timeout pattern is viable.
* **US-003:** As a developer, I want to generate a single valid conformer for cisplatin (SPL geometry) using Molassembler so that I can confirm the SCINE stereopermutator API is functional.
* **US-004:** As a developer, I want to record the exact Molassembler API surface used (module path, class names, method signatures) so that `MiniPRD_MolassemblerAdapter` has a concrete implementation target.

---

## 3. Implementation Plan (Task List)

- [x] Task 1: Run `uv add "scine-molassembler>=2.0.0"` to install SCINE into the project.
- [x] Task 2: Create `tests/spike_molassembler.py` (standalone script, not a unit test, not imported by main library).
- [x] Task 3: In `spike_molassembler.py`, import `scine.molassembler` (or the correct package path — verify against installed package) and print the version. Record the correct import path.
- [x] Task 4: Write a module-level worker function `_test_worker(args)` that does minimal Molassembler work. Call `pickle.dumps(_test_worker)` — record whether it succeeds or raises `PicklingError`.
- [x] Task 5: If picklable, submit `_test_worker` to `ProcessPoolExecutor(max_workers=1)` and call `fut.result(timeout=5)` — confirm the process isolation works.
- [x] Task 6: Write a minimal cisplatin conformer generation call: metal=Pt, geometry=SPL, 2×Cl + 2×NH₃ ligands. Log the output (XYZ block or error traceback). **[Value Error Fixed]**
- [x] Task 7: Save results to `.agents/memory/molassembler_spike_results.md` with: correct import path, picklability status, cisplatin conformer status, exact API method signatures used.

---

## 4. The Negative Space (Constraints)

* **DO NOT** write production code in this spike — all code lives in `tests/spike_molassembler.py` and is never imported by `src/oinsmiles/`.
* **DO NOT** assume `scine.molassembler` is the correct import path — verify against the installed package structure.
* **DO NOT** proceed to `MiniPRD_MolassemblerAdapter` until spike confirms picklability. If picklability fails, raise a new issue and propose `multiprocessing.Queue`-based alternative before proceeding.
* **DO NOT** commit `tests/spike_molassembler.py` as a permanent test — it is a temporary investigation artifact.

---

## 5. Integration Tests & Verification

* **Test 1 (Deterministic):** `import scine.molassembler` (or verified correct path) → no `ImportError`. **(PASS)**
* **Test 2 (Deterministic):** `pickle.dumps(_test_worker)` → no `PicklingError`. If this fails, document the exact error and halt. **(PASS)**
* **Test 3 (Deterministic):** `ProcessPoolExecutor(max_workers=1).submit(_test_worker, {}).result(timeout=5)` → returns without hanging. **(PASS)**
* **Test 4 (Novel):** Cisplatin DG conformer → Expected Output: [Candidate Artifact routing protocol triggered — XYZ block saved to `tests/candidate_outputs/spike_cisplatin.xyz` for human review. Do NOT assert specific atom coordinates.] **(PASS)**
