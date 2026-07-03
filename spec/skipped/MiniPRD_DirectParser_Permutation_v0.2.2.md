# MiniPRD: Direct Parser — Permutation Selection (v0.2.2 Blocker #3)

**Hypergraph Node ID:** `module_direct_parser_permutation_v022`
**Parent Node:** `system_direct_parser_v022`
**Parent SuperPRD:** `spec/compiled/SuperPRD_DirectParser_v0.2.2.md`
**Audit Reference:** `spec/audit/DirectParser_IntegrationAudit_20260506.md` — "Blocker #3: No Permutation/Isomerism Selection"; "Example 3: Legacy Adapter Comparison"
**Execution Order:** 2 of 5
**Blocked By:** MiniPRD #1 (FragmentMapping)
**Estimated Effort:** 6–8 hours (largest of the chain)
**Priority:** P0 (Block Ship)

---

## 1. The Confidence Mandate

**Confidence Score:** 7/10 (lowest of the chain; mandatory write-up gate raises it to ≥8 before code begins)

**Rationale:**
- ✅ Legacy `_pick_masm_permutation()` exists and is battle-tested in production.
- ✅ Audit Example 3 catalogs the symmetric-fixture cases (Cisplatin/Transplatin, fac/mer-Ir(ppy)₃) that drive permutation selection.
- ⚠️ Nobody has yet read `_pick_masm_permutation()` end-to-end with intent to port. Risk: state-coupling to `MolassemblerAdapter` internals.
- ⚠️ "Semantically equivalent port" is harder to verify than "byte-identical port" — see Test 1.

**Gating Deliverable (D4):** Before any code change, the implementer attaches a 1-page write-up to this MiniPRD as `spec/active/notes/permutation_legacy_readthrough.md` covering:
1. Call graph of `_pick_masm_permutation()` (every helper it calls).
2. State dependencies (which `MolassemblerAdapter` instance attributes it reads/writes).
3. Pure-vs-stateful classification of each branch.
4. Recommended port strategy: copy-as-pure-function vs. extract-shared-helper vs. redesign.

**Confidence is permitted to rise to 8/10 once the write-up lands.** If the write-up reveals deep state coupling (e.g., the function depends on a cached `Mol` from earlier in the legacy pipeline), MiniPRD #2 is paused and Architect is notified — possibly converting the MiniPRD into an extract-shared-module redesign.

---

## 2. Atomic User Stories

- **US-001:** As `parse_oin_direct()`, I want to select the correct Molassembler permutation index for symmetric metal complexes, so cis/trans, fac/mer, and other isomerism is preserved through the OIN→XYZ pipeline.
- **US-002:** As a developer maintaining the direct backend, I want the permutation logic in its own module (`src/oinsmiles/generation/permutation.py` or a clearly-bounded section of `engine.py`), so it can be audited independently of the rest of the parser.
- **US-003:** As a CI fixture suite, I want unit-level evidence that the direct backend's permutation selection produces the same *index value* as the legacy adapter on at least one canonical octahedral and one canonical TBP case.
- **US-004:** As an architect reviewing the port, I want a written read-through of the legacy implementation attached to this MiniPRD, so the port decision is evidence-based rather than trust-based.

---

## 3. Implementation Plan (Task List)

- [ ] **Task 1 (D4 gate — write-up):** Create `spec/active/notes/permutation_legacy_readthrough.md`. Read `MolassemblerAdapter._pick_masm_permutation()` end-to-end. Document call graph, state deps, pure-vs-stateful per branch. Conclude with port strategy recommendation. **Code may not begin until this lands.**
- [ ] **Task 2 (port strategy decision):** Based on Task 1's write-up, choose: (a) copy as pure free function, (b) extract a shared module both backends can import, or (c) redesign. If (c), pause and notify Architect. (a) and (b) proceed.
- [ ] **Task 3 (module skeleton):** Create `src/oinsmiles/generation/permutation.py` with `select_permutation_index(geometry_code, vertex_data, frag_to_atom_map) -> int`. Add `__all__ = ["select_permutation_index"]` for clean imports. (If write-up shows the function is short enough, fold into `engine.py` instead — A4 allows either.)
- [ ] **Task 4 (port — symmetric cases):** Implement permutation selection for: SPL (square-planar; Cisplatin/Transplatin), TBP (trigonal bipyramidal), OCT (octahedral; fac/mer-Ir(ppy)₃). Mirror legacy logic semantically. Add inline comments referencing the relevant Task-1 write-up section.
- [ ] **Task 5 (port — remaining geometries):** Implement for the remaining 6 geometry codes in `SUPPORTED_GEOMETRY_CODES`. For codes the legacy adapter handles via templates rather than permutation selection, document why permutation index is N/A and return a sentinel (e.g., `-1`) handled by the caller.
- [ ] **Task 6 (unsupported-geometry guard, R9):** If `geometry_code` is not in `SUPPORTED_GEOMETRY_CODES`, raise `NotImplementedError(f"Geometry {geometry_code} not supported by direct backend; try OIN3DGenerator(backend='legacy')")`.
- [ ] **Task 7 (parse_oin_direct integration):** Wire `select_permutation_index(...)` into `parse_oin_direct` after `_translate_vertex_indices` and before `construct_molassembler_mol`. Pass the result into Molassembler's permutation API.
- [ ] **Task 8 (unit test — Cisplatin/Transplatin parity):** `tests/unit/test_permutation.py::test_cisplatin_index_matches_legacy`. Construct the same input as the legacy adapter would receive; assert `select_permutation_index(...)` returns the same integer as `MolassemblerAdapter._pick_masm_permutation()` on Cisplatin and Transplatin. (Calls into `MolassemblerAdapter` for ground truth — read-only, no mutation.)
- [ ] **Task 9 (unit test — fac/mer-Ir(ppy)₃):** Same shape as Task 8 for octahedral fac/mer cases. Assert different integers for fac vs. mer.
- [ ] **Task 10 (cross-backend RMSD smoke):** Run `verify_roundtrip.py --backend direct` on the symmetric fixtures from Tasks 8–9. Mean RMSD < 1.0Å vs. input XYZ. (Full fixture suite verification is MiniPRD #5.)

---

## 4. The Negative Space (Constraints)

- **DO NOT** modify `MolassemblerAdapter._pick_masm_permutation()` or any other `molassembler_adapter.py` function. The legacy implementation is byte-stable per P6.4.
- **DO NOT** mutate `MolassemblerAdapter` instance state from the new module. If the port needs Adapter-internal state, that's a redesign signal — pause and notify Architect.
- **DO NOT** import `MolassemblerAdapter` anywhere in `permutation.py` (production code). Tests *may* import it for ground-truth comparison.
- **DO NOT** introduce a registry-of-permutation-strategies abstraction. One geometry-code dispatch is sufficient (YAGNI).
- **DO NOT** silently fall back from "permutation lookup failed" to "use default permutation 0." Failures must propagate as `NotImplementedError` or `ValueError` with the geometry code in the message.
- **DO NOT** start coding before the Task 1 write-up lands.
- **DO NOT** change `parse_oin_direct`'s signature in this MiniPRD; integration is a one-liner insertion.

---

## 5. Integration Tests & Verification

### Test 1 (Deterministic — Cisplatin index parity):
- **Input:** Cisplatin OIN string, with `frag_to_atom` from MiniPRD #1.
- **Expected Output:** `select_permutation_index(...)` returns the same `int` as legacy `MolassemblerAdapter._pick_masm_permutation()` on the same input. **This is the key port-fidelity assertion.**

### Test 2 (Deterministic — Transplatin):
- **Input:** Transplatin OIN string.
- **Expected Output:** Different `int` than Cisplatin (cis ≠ trans). Matches legacy.

### Test 3 (Deterministic — fac/mer-Ir(ppy)₃):
- **Input:** fac- and mer-Ir(ppy)₃ OIN strings.
- **Expected Output:** Different `int` per isomer. Matches legacy. Confirms octahedral discrimination.

### Test 4 (Cross-backend RMSD — round-trip):
- **Input:** Cisplatin, Transplatin, fac-Ir(ppy)₃, mer-Ir(ppy)₃ XYZ files.
- **Expected Output:** `verify_roundtrip.py --backend direct` mean RMSD < 1.0Å vs. input XYZ for each fixture.

### Test 5 (Unsupported geometry):
- **Input:** Synthetic OIN string with geometry code `XYZ` (not in `SUPPORTED_GEOMETRY_CODES`).
- **Expected Output:** `NotImplementedError` raised; message contains `"XYZ"` and the recommendation to try `backend='legacy'`.

### Candidate Artifact Routing:
- If permutation port produces a *different* index value than legacy on a symmetric fixture but the resulting mean RMSD is still < 1.0Å, route to `tests/fixtures/oin/_exclusions.yml` *only* with Architect approval (this is a "valid alternative permutation" case — easy to mistake for correctness; hard to ship without scrutiny).
- Any case where direct produces wrong isomerism (cis becomes trans, fac becomes mer) is a **hard block** — investigate root cause; do not exclude.
