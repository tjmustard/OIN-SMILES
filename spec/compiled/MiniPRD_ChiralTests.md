# MiniPRD: Chiral Test Suite

**Hypergraph Node ID:** `atom_ChiralTests` (new — tests module)
**Parent Node:** `sys_oin_smiles`

---

## 1. The Confidence Mandate

Before generating any plans or writing code, analyze this document and output a Confidence Score (1-10). If the score is below 9, list strictly the clarifying questions needed to reach 10.

**Prerequisites (all must be complete before this MiniPRD begins):**
1. `MiniPRD_ChiralEncoding` complete — `CIPAssigner`, `ChiralityRecoveryUtility`, `parse_inline_string()` fix all implemented and unit-tested.
2. `MiniPRD_MolassemblerAdapter` complete — all 6 existing OIN stability tests passing.
3. Candidate fixture SMILES have been generated (routed to `tests/candidate_outputs/`) and **human-reviewed** before promotion to `tests/fixtures/`.

**Candidate Artifact Protocol:** All new P-chiral, N-chiral, and axial-chiral SMILES fixture strings are novel outputs that must be saved to `tests/candidate_outputs/` for human review before any `assert` can be written against them in unit tests. The expected `@`/`@@` literals in Test 3 and Test 4 below come from the human-approved `tests/fixtures/` versions only.

---

## 2. Atomic User Stories

* **US-001:** As a QA process, I want BINAP to pass a structural stability test (pipeline does not crash, returns a non-None OIN string) without any `@`/`@@` P-atom assertion, so that the pipeline handles axial-chiral biaryl ligands gracefully.
* **US-002:** As a QA process, I want a P-chiral fixture (e.g., a chiral monodentate phosphine) to pass the RDKit CIP oracle test so that P stereocenter encoding correctness is independently verified (not just round-trip stable).
* **US-003:** As a QA process, I want an N-chiral fixture to pass the RDKit CIP oracle test so that N stereocenter encoding is also covered.
* **US-004:** As a QA process, I want an axial-chiral ligand (e.g., a BINAP-containing complex) where the SMILES output correctly encodes the axial chirality descriptor for those atoms, so that the encoded string is chemically correct even when we do not interpret axial chirality specially.
* **US-005:** As a QA process, I want all 6 existing regression complexes (cisplatin, transplatin, cis-PtCl₂(en), ferrocene, fac-/mer-Ir(ppy)₃) to still pass OIN stability after the full refactor, so that no regressions are introduced.

---

## 3. Implementation Plan (Task List)

- [ ] Task 1: Confirm all candidate fixture files in `tests/candidate_outputs/` have been human-reviewed. Promote approved fixtures to `tests/fixtures/` with their expected `@`/`@@` SMILES annotated in a companion `.expected` file or as a constant in the test file.
- [ ] Task 2: Write `tests/test_regression_stability.py` — run all 6 existing OIN stability checks post-refactor:
  - cisplatin, transplatin, cis-PtCl₂(en), ferrocene, fac-Ir(ppy)₃, mer-Ir(ppy)₃.
  - Assert: `encode(xyz_path)` == previously-approved OIN string (from `tests/fixtures/`).
- [ ] Task 3: Write `tests/test_binap_stability.py`:
  - `test_binap_does_not_crash()`: `assert XYZToSMILES().convert(BINAP_XYZ_PATH) is not None`.
  - **No** `@`/`@@` assertion on BINAP P atoms (BINAP chirality is axial, not P-centered).
- [ ] Task 4: Write `tests/test_chiral_p.py` with two test methods:
  - `test_p_stability()`: encode P-chiral fixture XYZ → OIN; decode → OIN'; assert `oin == oin'` (round-trip stable).
  - `test_p_cip_oracle()`: encode P-chiral fixture XYZ → OIN; call `extract_ligand_smiles(oin)` to strip metal/slots; call `Chem.MolFromSmiles(smiles)`; call `Chem.AssignStereochemistry(mol, cleanIt=True, force=True)`; assert `get_chiral_atom(mol, atomic_num=15).GetPropsAsDict().get('_CIPCode') == EXPECTED_DESCRIPTOR` (where `EXPECTED_DESCRIPTOR` is the human-verified R or S from `tests/fixtures/`).
- [ ] Task 5: Write `tests/test_chiral_n.py` — identical structure to `test_chiral_p.py` but using N-chiral fixture and `atomic_num=7`.
- [ ] Task 6: Write `tests/test_axial_chiral.py`:
  - For an axial-chiral ligand (e.g., BINAP or similar atropisomeric compound): encode XYZ → OIN.
  - Assert that the SMILES for the axial-chiral atoms contains the correct chirality descriptor characters (e.g., `/@` or `/` for axial chirality E/Z equivalents) as specified in the human-reviewed fixture.
  - This tests that the encoded SMILES is chemically correct for those atoms, not that OINInlineHandler interprets axial chirality.
- [ ] Task 7: Write a helper `extract_ligand_smiles(oin_string: str) -> str` in `tests/test_helpers.py` that strips the metal fragment and slot markers from an OIN string and returns the ligand SMILES for CIP oracle testing. (This is a test-only utility, not a production function.)
- [ ] Task 8: Run full test suite: `uv run python -m unittest discover tests`. All tests must pass (including existing regressions + new chiral tests).

---

## 4. The Negative Space (Constraints)

* **DO NOT** assert `@`/`@@` on BINAP P atoms — BINAP chirality is axial (atropisomeric), not P-centered. `test_binap_stability.py` asserts non-None only.
* **DO NOT** promote fixture SMILES from `tests/candidate_outputs/` to `tests/fixtures/` without human review and RDKit CIP oracle verification.
* **DO NOT** write CIP oracle tests using only round-trip (OIN-in == OIN-out) — the RDKit CIP oracle cross-validator (`test_p_cip_oracle`, `test_n_cip_oracle`) is mandatory for correctness verification.
* **DO NOT** read from `tests/candidate_outputs/` in CI — those are unverified outputs, not regression baselines.
* **DO NOT** write `extract_ligand_smiles()` as a production function in `src/oinsmiles/` — it is a test helper only.
* **DO NOT** hardcode expected `@`/`@@` SMILES literals in test files without annotating the source (which XYZ, what date, who reviewed).

---

## 5. Integration Tests & Verification

* **Test 1 (Deterministic):** All 6 existing OIN stability checks pass after MolassemblerAdapter refactor: cisplatin, transplatin, cis-PtCl₂(en), ferrocene, fac-/mer-Ir(ppy)₃.
* **Test 2 (Deterministic):** `test_binap_does_not_crash()` → `XYZToSMILES().convert(BINAP_XYZ_PATH) is not None`. No `@`/`@@` assertion.
* **Test 3 (Deterministic):** `test_p_cip_oracle()` → RDKit `_CIPCode` on P atom == human-verified expected descriptor (R or S). This is the independent correctness oracle that breaks the circular validation.
* **Test 4 (Deterministic):** `test_n_cip_oracle()` → RDKit `_CIPCode` on N atom == human-verified expected descriptor.
* **Test 5 (Novel):** Axial-chiral ligand SMILES correctness check → [Candidate Artifact routing triggered — annotated SMILES string with axial chirality descriptors saved to `tests/candidate_outputs/axial_chiral_encoded.smi` for human review before `test_axial_chiral.py` assertions can be written]
