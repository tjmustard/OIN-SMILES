# MiniPRD: Direct Parser — Eta-Bond Atom-Index Translation (v0.2.2 Blocker #4)

**Hypergraph Node ID:** `module_direct_parser_eta_bonds_v022`
**Parent Node:** `system_direct_parser_v022`
**Parent SuperPRD:** `spec/compiled/SuperPRD_DirectParser_v0.2.2.md`
**Audit Reference:** `spec/audit/DirectParser_IntegrationAudit_20260506.md` — "Blocker #4: Eta Bond Translation Broken"
**Execution Order:** 3 of 5
**Blocked By:** MiniPRD #1 (FragmentMapping)
**Estimated Effort:** 4–6 hours
**Priority:** P1 (Pre-merge)

---

## 1. The Confidence Mandate

**Confidence Score:** 8/10

**Rationale:**
- ✅ MiniPRD #1's `frag_to_atom` mapping is the missing input that makes eta translation tractable.
- ✅ `TEMPLATES` data in `oin_aligner.py` already encodes slot resolution for η⁵ rings; reusable.
- ✅ Existing `oin_parser.py:162` `IndexError` catch (per CLAUDE.md memory) suggests the failure mode is well-understood.
- ⚠️ Multi-eta cases (ansa-metallocenes like TiCat1/3/4) need careful index management — easy to off-by-one.
- ⚠️ "No bond-to-self" is a narrow regression test that won't catch all mis-translations; supplemented by RMSD round-trip.

---

## 2. Atomic User Stories

- **US-001:** As `parse_oin_direct()`, I want to translate η-bond vertex indices (in OIN's ring-rank space) to atom indices (in connected-SMILES space) using `frag_to_atom`, so ferrocene and TiCp₂Me₂ generate correct sandwich geometry.
- **US-002:** As a multi-eta complex (TiCat1/3/4), I want each η⁵ ring's atoms identified and bonded correctly to the metal, so ansa-metallocenes round-trip without losing ring connectivity.
- **US-003:** As a maintainer, I want a named regression test pinning the exact "bond to self" failure from the audit, so the bug cannot silently re-occur.

---

## 3. Implementation Plan (Task List)

- [ ] **Task 1 (helper function):** Add `_translate_eta_vertex_to_atoms(vertex_index: int, fragment_rank: int, frag_to_atom: Dict[int, List[int]], template_spec: Dict) -> List[int]` in `engine.py` (or `polydentate.py` if already created in MiniPRD #4 — coordinate at module-creation time).
- [ ] **Task 2 (single-eta translation):** For η⁵ Cp on ferrocene, given a fragment rank and a vertex index, return the atom indices of the 5 ring atoms. Use `frag_to_atom[fragment_rank]` to get the candidate atom set; use `TEMPLATES`/`TEMPLATE_SPECS` to identify which subset is the ring.
- [ ] **Task 3 (multi-eta translation):** Extend Task 2 to handle complexes with multiple η-ligands (TiCat1/3/4 ansa-metallocenes). Two η⁵ rings → two distinct slot resolutions, each producing 5 atom indices.
- [ ] **Task 4 (parse_oin_direct wiring):** In the eta-bond construction step of `parse_oin_direct`, replace `vertex_index` arguments to Molassembler bond-add calls with `atom_index` arguments produced by Task 1's helper.
- [ ] **Task 5 (bond-to-self prevention):** Add a precondition assertion: `assert source_atom_idx != target_atom_idx, f"Refusing bond-to-self at atom {source_atom_idx}"`. This is the named regression-test target (A9).
- [ ] **Task 6 (regression test, A9):** Create `tests/unit/test_eta_bonds.py::test_no_bond_to_self_after_eta_translation`. Constructs an OIN string for ferrocene; calls `parse_oin_direct`; asserts no Molassembler bond has equal source and target atom indices.
- [ ] **Task 7 (round-trip — single-eta):** Add `tests/unit/test_eta_bonds.py::test_ferrocene_roundtrip`. `verify_roundtrip.py --backend direct` on `tests/fixtures/oin/ferrocene.json`. Mean RMSD < 1.0Å.
- [ ] **Task 8 (round-trip — multi-eta):** Add `test_ticat1_roundtrip`, `test_ticat3_roundtrip`, `test_ticat4_roundtrip`. Each runs through direct backend; asserts mean RMSD < 1.0Å.
- [ ] **Task 9 (multi-eta atom-distinctness test):** Add `test_ansa_metallocene_rings_disjoint`. Constructs a TiCat1 OIN string; asserts the two η⁵ rings have *disjoint* atom-index sets (no atom belongs to both rings). Catches the failure mode where multi-eta translation conflates two rings.

---

## 4. The Negative Space (Constraints)

- **DO NOT** modify `TEMPLATES` or `TEMPLATE_SPECS` in `oin_aligner.py` or anywhere else; read-only consumption.
- **DO NOT** modify `MolassemblerAdapter` or `_stitch_multi_eta_fragment` (legacy backend dependency; P6.4).
- **DO NOT** silently skip eta bonds when translation fails. Raise `OINFormatError(f"Failed to translate eta vertex {v} for fragment {f}: {detail}")`.
- **DO NOT** assume η⁵ — handle the η-hapticity from the OIN string. Other haptic modes (η³, η⁶, η⁸) must produce a clear error if encountered (`NotImplementedError`) until a future MiniPRD adds them.
- **DO NOT** add a "best-guess" fallback that picks the first N atoms in a fragment if template-spec lookup fails. Hard-fail with a diagnostic message instead.
- **DO NOT** modify `oin_parser.py:162` IndexError catch (already adjusted per CLAUDE.md memory; that's separate scope).

---

## 5. Integration Tests & Verification

### Test 1 (Deterministic — bond-to-self regression, A9):
- **Input:** Ferrocene OIN string.
- **Expected Output:** `parse_oin_direct(...)` completes; no Molassembler bond has `source == target`.

### Test 2 (Deterministic — ferrocene round-trip):
- **Input:** `tests/fixtures/oin/ferrocene.json`.
- **Expected Output:** `verify_roundtrip.py --backend direct` returns mean RMSD < 1.0Å vs. input XYZ.

### Test 3 (Deterministic — TiCat1/3/4 round-trip):
- **Input:** `tests/fixtures/oin/ticat{1,3,4}.json`.
- **Expected Output:** All three pass `verify_roundtrip.py --backend direct` with mean RMSD < 1.0Å.

### Test 4 (Deterministic — disjoint multi-eta rings):
- **Input:** TiCat1 OIN string.
- **Expected Output:** The atom-index sets returned by `_translate_eta_vertex_to_atoms` for ring A and ring B are disjoint (no shared atoms).

### Test 5 (Negative — unsupported hapticity):
- **Input:** Synthetic OIN string declaring η⁶-arene.
- **Expected Output:** `NotImplementedError` raised; message names the hapticity and points to `backend='legacy'`.

### Candidate Artifact Routing:
- A new η-ligand fixture introduced post-v0.2.2 that fails Test 1's bond-to-self check is a **hard block** — investigate; do not exclude.
- Mean RMSD between 1.0Å and 1.5Å on a TiCat fixture routes to `_exclusions.yml` *with documented justification* (e.g., "ring-pucker variance under DG ensemble"); requires Architect approval.
