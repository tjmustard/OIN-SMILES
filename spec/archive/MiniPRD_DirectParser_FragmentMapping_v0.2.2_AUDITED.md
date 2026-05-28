# MiniPRD: Direct Parser — Fragment Mapping (v0.2.2 Blocker #1)

**Hypergraph Node ID:** `module_direct_parser_fragment_mapping_v022`
**Parent Node:** `system_direct_parser_v022`
**Parent SuperPRD:** `spec/compiled/SuperPRD_DirectParser_v0.2.2.md`
**Audit Reference:** `spec/audit/DirectParser_IntegrationAudit_20260506.md` — "Blocker #1: Fragment Rank ↔ Atom Index Mapping"
**Execution Order:** 1 of 5 (foundation; all subsequent MiniPRDs depend on this)
**Estimated Effort:** 2–3 hours
**Priority:** P0 (Block Ship)

---

## 1. The Confidence Mandate

**Confidence Score:** 9/10

**Rationale:**
- ✅ Audit document Example 1 ("Mapping Problem Illustrated") is the primary reference; the failure mode is fully understood.
- ✅ Metal-first canonical-form invariant codified in project memory; no ambiguity about `fragments[0]`.
- ✅ Pre-MiniPRD-#1 verification spike (per D3) gates start: ~30-line script confirms `fragments[0]` is the metal across every fixture in `tests/fixtures/oin/*.json` before any code change lands.
- ⚠️ The connected-SMILES builder must produce a deterministic atom ordering across runs; one helper test pins this.

**Gating Deliverable (D3):** Run `tools/verify_metal_first.py` (created in Task 1) over all fixtures *before* committing Task 3. If any fixture returns a non-metal at `fragments[0]`, halt and escalate to Architect — that fixture is malformed per project invariant.

---

## 2. Atomic User Stories

- **US-001:** As `_extract_oin_constraints()`, I want to return a `Dict[int, List[int]]` mapping fragment-rank → list of atom-indices in the connected SMILES, so downstream stages (eta-bond translation, polydentate connectivity) can convert OIN's vertex-rank coordinates to RDKit/Molassembler atom indices.
- **US-002:** As callers of `_extract_oin_constraints()`, I want the function to be underscore-prefixed (private), so the v0.2.1 `extract_oin_constraints` re-export from `engine.py` becomes an internal-only helper without a public-API breaking-change announcement.
- **US-003:** As a developer running the test suite, I want a script that enumerates every call site of `_extract_oin_constraints` (in src/ and tests/) and produces a count, so the "17 call sites" claim is replaced with live data.
- **US-004:** As a project invariant guardian, I want a verification spike that asserts `fragments[0]` is a transition metal across every existing fixture, so we run with documented evidence rather than implicit trust.

---

## 3. Implementation Plan (Task List)

- [ ] **Task 1 (verification spike, D3 gate):** Create `tools/verify_metal_first.py`. For every `tests/fixtures/oin/*.json`, parse the OIN string, split into fragments, assert fragment[0] is a transition-metal symbol. Print per-fixture pass/fail. Run before any code change. Time-box: 30 minutes.
- [ ] **Task 2 (call-site audit tool, A8):** Create `tools/audit_extract_calls.py`. `grep`s for `extract_oin_constraints` across `src/` and `tests/`, prints a one-per-line list of `path:line: snippet`. Will be re-run after rename to confirm zero leftover unprefixed references.
- [ ] **Task 3 (rename to private):** In `src/oinsmiles/generation/engine.py`, rename `extract_oin_constraints` → `_extract_oin_constraints`. Confirm it is *not* in `__all__` of either `engine.py` or `__init__.py` (per project memory it never was; this is bookkeeping confirmation, not a public-API change).
- [ ] **Task 4 (signature change):** Change return type from `Tuple[str, Dict[int, Dict]]` to `Tuple[str, Dict[int, Dict], Dict[int, List[int]]]`. The third element is the `fragment_to_atom_mapping`.
- [ ] **Task 5 (mapping construction):** Inside `_extract_oin_constraints`, after the connected-SMILES builder produces the joined string, walk the connected SMILES and record which atom indices originated from which fragment rank. Use the regex preprocessor's existing fragment-split state as ground truth.
- [ ] **Task 6 (caller update — production):** Update the single production caller (`parse_oin_direct` in `engine.py`) to unpack the new 3-tuple. Use a `_, _, frag_to_atom = _extract_oin_constraints(s)` style for forward compatibility, but capture the new mapping into a named local for use in subsequent MiniPRDs.
- [ ] **Task 7 (caller update — tests):** Run `tools/audit_extract_calls.py`; update every test call site to either unpack the 3-tuple or use `_extract_oin_constraints(s)[:2]` if the third element is unused. Re-run the tool to verify zero unprefixed references remain.
- [ ] **Task 8 (deterministic-ordering test):** Add `tests/unit/test_fragment_mapping.py::test_atom_ordering_deterministic`. Calls `_extract_oin_constraints` 10× on the same Cisplatin OIN string; asserts `frag_to_atom` is byte-identical across all 10 calls.
- [ ] **Task 9 (round-trip smoke):** Add `tests/unit/test_fragment_mapping.py::test_cisplatin_mapping_correct`. Asserts `frag_to_atom[0] == [0]` (metal at fragment 0, atom 0) and `frag_to_atom[1]` and `frag_to_atom[2]` each contain exactly one Cl atom index. No "bond to self" — a separate canary in MiniPRD #5.
- [ ] **Task 10 (CHANGELOG):** Add an entry under "v0.2.2 — Unreleased" noting that `extract_oin_constraints` has been privatized to `_extract_oin_constraints` and gained a third return element. Mention it is a rename of an undocumented helper, not a public-API change.

---

## 4. The Negative Space (Constraints)

- **DO NOT** make `_extract_oin_constraints` public again or add it to `__all__`.
- **DO NOT** change the connected-SMILES builder's output beyond what's required to expose the mapping (zero functional change to the resulting SMILES).
- **DO NOT** modify `MolassemblerAdapter` or anything in `molassembler_adapter.py` (P6.4).
- **DO NOT** touch eta-bond logic, polydentate logic, or permutation logic in this MiniPRD; they are dedicated MiniPRDs #3, #4, #2 respectively.
- **DO NOT** add a runtime metal-detection branch — the metal-first invariant holds. (If Task 1 finds a counterexample, halt and escalate; do not patch around it.)
- **DO NOT** change `tokenize_unsanitized_smiles()` in this MiniPRD; the upstream API is stable.
- **DO NOT** add new dependencies.

---

## 5. Integration Tests & Verification

### Test 1 (Deterministic — Cisplatin):
- **Input:** OIN string for Cisplatin (`tests/fixtures/oin/cisplatin.json`).
- **Expected Output:** `_extract_oin_constraints(s)` returns 3-tuple. `frag_to_atom[0] == [0]` (Pt at atom 0). `frag_to_atom[1]` and `frag_to_atom[2]` each have one element. Atom indices are contiguous and cover the full atom count.

### Test 2 (Deterministic — ordering stability):
- **Input:** Run `_extract_oin_constraints` 10 times on `cisplatin.json`.
- **Expected Output:** Identical `frag_to_atom` dict on every call (no dict-ordering noise).

### Test 3 (Deterministic — polydentate sanity):
- **Input:** OIN string for `cis-PtCl2(en)` (en = ethylenediamine, bidentate).
- **Expected Output:** `frag_to_atom[1]` has ≥2 atom indices (the two N atoms of en plus their carbons; full polydentate handling is MiniPRD #4 territory but the *mapping* should already include all atoms).

### Test 4 (Verification spike — gating):
- **Input:** Every `tests/fixtures/oin/*.json` parsed by `tools/verify_metal_first.py`.
- **Expected Output:** All fixtures pass. **If any fixture fails, MiniPRD #1 is halted and the failure routed to Architect (per D3 gate).**

### Test 5 (Call-site audit):
- **Input:** `tools/audit_extract_calls.py` run after Task 7.
- **Expected Output:** Zero matches for unprefixed `extract_oin_constraints` (confirming rename is complete). Match list for `_extract_oin_constraints` contains the expected production caller and all tests.

### Candidate Artifact Routing:
- A novel test fixture introduced post-v0.2.2 that fails Task 1's metal-first check is routed to `tests/fixtures/oin/_exclusions.yml` with justification (per A12); does not block MiniPRD #1 close.
- Any deterministic-ordering test failure is a **hard block** — atom indices must be byte-stable.
