# MiniPRD: Direct Parser — Polydentate Ligand Connectivity (v0.2.2 Blocker #2)

**Hypergraph Node ID:** `module_direct_parser_polydentate_v022`
**Parent Node:** `system_direct_parser_v022`
**Parent SuperPRD:** `spec/compiled/SuperPRD_DirectParser_v0.2.2.md`
**Audit Reference:** `spec/audit/DirectParser_IntegrationAudit_20260506.md` — "Blocker #2: Polydentate Ligand Connectivity"
**Execution Order:** 4 of 5
**Blocked By:** MiniPRD #1 (FragmentMapping), MiniPRD #2 (Permutation)
**Estimated Effort:** 4–6 hours
**Priority:** P1 (Pre-merge)

---

## 1. The Confidence Mandate

**Confidence Score:** 7/10

**Rationale:**
- ✅ Existing chiral fixtures (PdCl₂-RR-BDPP, PdCl₂-R-BINAP, PdCl₂-RR-BDNN) provide ground truth for polydentate behavior.
- ✅ `vector_data` from `parse_inline_string()` already exposes binding-atom information; this MiniPRD wires it through the direct path.
- ⚠️ Bidentate vs. tridentate vs. tetradentate may need different code paths; one architectural decision to make during implementation.
- ⚠️ Per the SuperPRD risk register: "polydentate mapping requires architectural change to vector data" is a *Medium / High* risk item — escalation path is via Architect.

---

## 2. Atomic User Stories

- **US-001:** As `parse_oin_direct()`, I want to identify all binding atoms of a polydentate ligand and add separate metal-ligand bonds for each, so chelating ligands like en, BINAP, BDPP, BDNN produce correct connectivity.
- **US-002:** As a chiral diphosphine fixture (BDPP, BDNN), I want my P-stereocenters preserved through the direct backend, so chirality round-trips with mean RMSD < 1.0Å.
- **US-003:** As a maintainer, I want polydentate logic encapsulated in a clearly-bounded function (in `engine.py` or `polydentate.py` if size warrants), so future ligand types can be added without sprawling edits.

---

## 3. Implementation Plan (Task List)

- [ ] **Task 1 (binding-atom identification):** Add `_identify_binding_atoms(fragment_rank, vector_data, frag_to_atom) -> List[int]`. Returns the atom indices in the fragment that bind directly to the metal. Use `vector_data` (from `parse_inline_string`) as authoritative; cross-check against `frag_to_atom`.
- [ ] **Task 2 (multi-bond construction):** In `parse_oin_direct`, for each polydentate fragment, iterate over `_identify_binding_atoms(...)` and add one Molassembler metal-ligand bond per binding atom (instead of one bond per fragment).
- [ ] **Task 3 (denticity dispatch):** If `len(binding_atoms) == 1`: monodentate path (current behavior). If `== 2`: bidentate. If `>= 3`: polydentate. Behavior is uniform (each binding atom gets its own bond), but the dispatch doc-strings the cases for clarity.
- [ ] **Task 4 (chirality preservation check):** For BDPP/BDNN fixtures, ensure the @/@@ chiral tags on the binding P/N atoms (preserved through MiniPRD #1 and #2) survive into the Molassembler molecule. If they don't, the issue lies in `tokenize_unsanitized_smiles` or earlier — file a follow-up issue; do not patch around in this MiniPRD.
- [ ] **Task 5 (en — Pt(en)Cl₂ round-trip):** `tests/unit/test_polydentate.py::test_pt_en_cl2_roundtrip`. `verify_roundtrip.py --backend direct` on `tests/fixtures/oin/cis-PtCl2-en.json`. Mean RMSD < 1.0Å.
- [ ] **Task 6 (BINAP — atropisomeric stability):** `test_binap_roundtrip`. PdCl₂-R-BINAP fixture; mean RMSD < 1.0Å.
- [ ] **Task 7 (BDPP — P-stereocenter):** `test_bdpp_p_chirality_roundtrip`. PdCl₂-RR-BDPP fixture; assert R/R configuration is preserved (CIP code matches input); mean RMSD < 1.0Å.
- [ ] **Task 8 (BDNN — N-stereocenter):** `test_bdnn_n_chirality_roundtrip`. PdCl₂-RR-BDNN fixture; analogous to Task 7.
- [ ] **Task 9 (binding-atom-count test):** `test_binding_atoms_count_correct`. For en, asserts `len(_identify_binding_atoms(en_fragment_rank, ...)) == 2`. For BINAP/BDPP, also `== 2`. For tridentate (if any future fixture), `== 3`.

---

## 4. The Negative Space (Constraints)

- **DO NOT** modify `parse_inline_string()` or its return contract. `vector_data` is consumed read-only.
- **DO NOT** modify the chirality module (`core/chirality.py`); polydentate handling does not change CIP assignment.
- **DO NOT** add a special-case branch per ligand type (en-specific, BINAP-specific code paths). Generic denticity dispatch via Task 3 is the only allowed shape.
- **DO NOT** change the metal-ligand bond order in Molassembler — denticity is encoded by *number of bonds*, not by bond-order multiplication.
- **DO NOT** silently skip a binding atom if its index isn't in `frag_to_atom`. Hard-fail with `OINFormatError(f"Polydentate binding atom {idx} not in fragment {rank} mapping; check vector_data")`.

---

## 5. Integration Tests & Verification

### Test 1 (Deterministic — Pt(en)Cl₂):
- **Input:** `tests/fixtures/oin/cis-PtCl2-en.json`.
- **Expected Output:** `parse_oin_direct(...)` produces a Molassembler molecule with 2 Pt–N bonds (one per N of en) and 2 Pt–Cl bonds. `verify_roundtrip.py --backend direct` mean RMSD < 1.0Å.

### Test 2 (Deterministic — BINAP):
- **Input:** PdCl₂-R-BINAP fixture.
- **Expected Output:** 2 Pd–P bonds (one per P of BINAP). Mean RMSD < 1.0Å.

### Test 3 (Deterministic — BDPP P-chirality):
- **Input:** PdCl₂-RR-BDPP fixture.
- **Expected Output:** 2 Pd–P bonds; round-tripped XYZ produces SMILES with @/@@ tags matching input (R/R configuration). Mean RMSD < 1.0Å.

### Test 4 (Deterministic — BDNN N-chirality):
- **Input:** PdCl₂-RR-BDNN fixture.
- **Expected Output:** 2 Pd–N bonds; R/R configuration preserved. Mean RMSD < 1.0Å.

### Test 5 (Deterministic — denticity correctness):
- **Input:** Each fixture above.
- **Expected Output:** `_identify_binding_atoms(...)` returns the expected count (2 for all four fixtures).

### Test 6 (Negative — missing binding atom):
- **Input:** Synthetic OIN string where `vector_data` references a binding atom not in `frag_to_atom`.
- **Expected Output:** `OINFormatError` raised with the missing index in the message.

### Candidate Artifact Routing:
- If P-stereocenter chirality is *lost* on round-trip but RMSD passes (R becomes S geometrically equivalent under symmetry), this is a **hard block** — chirality loss is a science-grade bug, not a tolerance issue. Investigate `CIPAssigner.assign_all()` interaction with the direct path before excluding.
- If RMSD is between 1.0Å and 1.5Å on a chiral diphosphine, route to `_exclusions.yml` *only* if chirality is preserved and the divergence is conformational (e.g., backbone-rotation variance); requires Architect approval.
