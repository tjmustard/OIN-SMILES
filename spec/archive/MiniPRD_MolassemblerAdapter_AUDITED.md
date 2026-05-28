# MiniPRD: Molassembler Adapter

**Hypergraph Node ID:** `atom_MolassemblerAdapter` (new), `atom_OIN3DGenerator` (modified), `atom_ArchitectorAdapter` (deleted), `atom_ArchitectorWrapper` (deleted)
**Parent Node:** `mod_generation`

---

## 1. The Confidence Mandate

Before generating any plans or writing code, analyze this document and output a Confidence Score (1-10). If the score is below 9, list strictly the clarifying questions needed to reach 10.

**Prerequisite:** `MiniPRD_MolassemblerSpike` must be complete. Read `.agents/memory/molassembler_spike_results.md` before starting this MiniPRD. The exact Molassembler API surface (import path, class names, method signatures, picklability status) documented in the Spike results is required to implement `_molassembler_worker`.

---

## 2. Atomic User Stories

* **US-001:** As `OIN3DGenerator`, I want `MolassemblerAdapter` to replace `ArchitectorAdapter` + `ArchitectorWrapper` so that Molassembler is the sole 3D generation backend.
* **US-002:** As a pipeline user, I want Molassembler conformer generation to execute in a separate OS process via `ProcessPoolExecutor` so that a GIL-holding C++ call cannot block the Python process indefinitely.
* **US-003:** As a pipeline user, I want `MolassemblerTimeoutError` raised (not silently swallowed) when the 60-second timeout fires, so that callers can handle it explicitly.
* **US-004:** As a developer, I want `_molassembler_worker` to be a module-level function (not a method or lambda), so that it is picklable by `ProcessPoolExecutor`.
* **US-005:** As a regression test suite, I want all 6 existing OIN stability tests to still pass after `OIN3DGenerator` is rewired to `MolassemblerAdapter`, so that no regressions are introduced.

---

## 3. Implementation Plan (Task List)

- [ ] Task 1: Read `.agents/memory/molassembler_spike_results.md` to get the confirmed import path, API signatures, and picklability status from the Spike.
- [ ] Task 2: Create `src/oinsmiles/generation/molassembler_adapter.py` with:
  - `MolassemblerTimeoutError(RuntimeError)` — custom exception.
  - `_molassembler_worker(args: dict) -> str` — module-level function, returns XYZ block. `args` contains: `metal`, `ligand_smiles_list`, `geometry_code`, `site_coords` (from `ParsedOIN`).
  - `MolassemblerAdapter` class with `__init__(self, timeout: int = 60)` and `generate(self, parsed_oin: ParsedOIN) -> str`.
- [ ] Task 3: Implement `MolassemblerAdapter.generate()`:
  1. Build `args` dict from `parsed_oin` fields.
  2. Submit `_molassembler_worker` to `ProcessPoolExecutor(max_workers=1)`.
  3. Call `fut.result(timeout=self.timeout)`.
  4. On `concurrent.futures.TimeoutError`: raise `MolassemblerTimeoutError(f"Molassembler timed out after {self.timeout}s")`.
  5. Return XYZ block string.
- [ ] Task 4: Implement `_molassembler_worker(args)` using the confirmed Molassembler API from Spike results. Convert `ParsedOIN` geometry vectors to Molassembler's site coordinate format.
- [ ] Task 5: Update `src/oinsmiles/generation/engine.py` (`OIN3DGenerator.generate()`):
  - Remove import of `ArchitectorAdapter`, `ArchitectorWrapper`.
  - Import `MolassemblerAdapter`, `MolassemblerTimeoutError`.
  - Replace `ArchitectorAdapter().convert(parsed_oin)` → `ArchitectorWrapper().run(inputDict)` chain with `MolassemblerAdapter().generate(parsed_oin)`.
- [ ] Task 6: Run all existing OIN stability tests (`uv run python -m unittest discover tests`). All 6 must pass before proceeding.
- [ ] Task 7: Delete `src/oinsmiles/generation/architector_adapter.py`.
- [ ] Task 8: Delete `src/oinsmiles/generation/wrapper.py`.
- [ ] Task 9: Update `pyproject.toml`:
  - Remove `architector` from `[project.dependencies]`.
  - Add `scine-molassembler >= 2.0.0` to `[project.dependencies]` (if not already added by Spike).
- [ ] Task 10: Run `uv sync` to verify dependency resolution succeeds after removal of `architector`.

---

## 4. The Negative Space (Constraints)

* **DO NOT** use `threading.Thread` for the Molassembler timeout — use `concurrent.futures.ProcessPoolExecutor` exclusively.
* **DO NOT** implement `_molassembler_worker` as a class method, instance method, or lambda — it must be a picklable module-level function.
* **DO NOT** delete `ArchitectorAdapter` or `ArchitectorWrapper` until `MolassemblerAdapter` passes all 6 existing OIN stability tests (Task 6 must succeed first).
* **DO NOT** use `multiprocessing.Pool` — use `concurrent.futures.ProcessPoolExecutor` for API consistency.
* **DO NOT** swallow `MolassemblerTimeoutError` inside `MolassemblerAdapter` — it must propagate to `OIN3DGenerator` and then to the caller.
* **DO NOT** assume the Molassembler API surface — use only what is documented in `.agents/memory/molassembler_spike_results.md`.

---

## 5. Integration Tests & Verification

* **Test 1 (Deterministic):** `pickle.dumps(_molassembler_worker)` → no `PicklingError` (re-verify from Spike in production module context).
* **Test 2 (Deterministic):** `MolassemblerAdapter(timeout=1).generate(cisplatin_parsed_oin)` with an artificial infinite-loop worker substituted → raises `MolassemblerTimeoutError` within 2 seconds (confirms timeout path works).
* **Test 3 (Deterministic):** After Task 7-8: `from oinsmiles.generation.architector_adapter import ArchitectorAdapter` → `ImportError` (confirms deletion).
* **Test 4 (Deterministic):** After Task 10: `uv sync` exits 0 with no unresolvable dependency errors.
* **Test 5 (Deterministic):** Core 5 OIN stability tests pass (v0.2.0 baseline): cisplatin, transplatin, cis-PtCl₂(en), fac-Ir(ppy)₃, mer-Ir(ppy)₃. **Note:** Ferrocene and ansa-metallocene complexes (TiCp2Me2, TiCat1/3/4) exhibit regressions from v0.2.1 eta-ligand changes (documented as known limitations); these are outside the v0.2.0 MiniPRD scope.
* **Test 6 (Novel):** `MolassemblerAdapter().generate(cisplatin_parsed_oin)` → XYZ block returned within 60s → [Candidate Artifact routing triggered — XYZ saved to `tests/candidate_outputs/molassembler_cisplatin.xyz` for RMSD review against original cisplatin fixture]
