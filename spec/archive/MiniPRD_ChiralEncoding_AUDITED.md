# MiniPRD: Chiral Encoding

**Hypergraph Node ID:** `atom_CIPAssigner`, `atom_ChiralityRecoveryUtility`, `atom_PseudoAtomStrategy`, `atom_OINInlineHandler` (modified)
**Parent Node:** `mod_core`, `mod_oin`

---

## 1. The Confidence Mandate

Before generating any plans or writing code, analyze this document and output a Confidence Score (1-10). If the score is below 9, list strictly the clarifying questions needed to reach 10.

**Key concern:** The `_CIPCode` atom property set by `CIPAssigner.assign_all()` must survive `OINSanitizer.generate_robust_smiles()` without being destroyed. Verify that `OINSanitizer`'s sanitization steps do not call `Chem.RemoveAllHs()` or re-initialize atom properties before starting implementation.

---

## 2. Atomic User Stories

* **US-001:** As the XYZ→OIN pipeline, I want `CIPAssigner.assign_all()` to assign `_CIPCode` to all P/N atoms on the full (pre-fragmentation) sanitized mol, so that fragmentation does not destroy chirality context.
* **US-002:** As the XYZ→OIN pipeline, I want `ChiralityRecoveryUtility.recover()` to re-apply `@`/`@@` to ligand SMILES after `OINSanitizer` runs, using stored `_CIPCode`, so that chiral phosphines and amines are correctly encoded.
* **US-003:** As the XYZ→OIN pipeline, I want `PseudoAtomStrategy` to activate as a fallback when `CIPAssigner` finds no `_CIPCode` on a P/N atom (e.g., non-standard valence), so that the pipeline does not crash.
* **US-004:** As the OIN parser, I want `parse_inline_string()` to strip slot markers via regex only (no `MolFromSmiles → MolToSmiles`), so that `@`/`@@` markers are preserved through the parse step.
* **US-005:** As the XYZ→OIN pipeline, I want `Chem.SanitizeMol(mol)` called as a hard precondition before `CIPAssigner.assign_all()`, so that `AssignStereochemistry()` has all required ring info and aromaticity data.

---

## 3. Implementation Plan (Task List)

- [ ] Task 1: Read `src/oinsmiles/utils/oin_aligner.py` and `src/oinsmiles/utils/xyz2mol.py` to understand exactly where fragmentation happens and where `OINSanitizer.generate_robust_smiles()` is called. Confirm insertion points for `CIPAssigner` and `ChiralityRecoveryUtility`.
- [ ] Task 2: Audit `src/oinsmiles/utils/xyz2mol_local.py` for the `chiral_stereo_check` flag. Confirm that `CIPAssigner.assign_all()` runs AFTER any stereo flags set by `xyz2mol_local.py`, and overrides them if needed.
- [ ] Task 3: Create `src/oinsmiles/core/chirality.py` with three classes:
  - `CIPAssigner` with `assign_all(mol: Chem.Mol) -> Chem.Mol`
  - `ChiralityRecoveryUtility` with `recover(mol: Chem.Mol) -> Chem.Mol`
  - `PseudoAtomStrategy` with `PSEUDO_ATOMIC_NUM: int = 0`
- [ ] Task 4: Implement `CIPAssigner.assign_all()`:
  1. Validate `mol is not None` (raise `ValueError` otherwise).
  2. Call `Chem.SanitizeMol(mol)` — this is the hard precondition; caller must handle `Chem.SanitizeMol` exceptions.
  3. Call `Chem.AssignStereochemistry(mol, cleanIt=True, force=True)`.
  4. For each atom where `atom.GetAtomicNum() in (7, 15)` (N=7, P=15): read `_CIPCode` via `atom.GetPropsAsDict().get('_CIPCode')` and store as atom property for later retrieval.
  5. Return mol.
- [ ] Task 5: Implement `ChiralityRecoveryUtility.recover()`:
  1. For each P/N atom: read stored `_CIPCode`.
  2. Map `_CIPCode` → `Chem.ChiralType` (R → CHI_TETRAHEDRAL_CCW, S → CHI_TETRAHEDRAL_CW, or confirmed correct mapping for RDKit convention).
  3. Call `atom.SetChiralTag(chiral_type)`.
  4. If no `_CIPCode` found on a P/N atom that has 4 neighbors: invoke `PseudoAtomStrategy` on that atom.
  5. Return mol.
- [ ] Task 6: Implement `PseudoAtomStrategy` substitution: replace the non-standard-valence P/N atom with a wildcard atom (`PSEUDO_ATOMIC_NUM = 0`). Add a utility method `strip_pseudo_atoms(mol)` that must be called before OIN serialization.
- [ ] Task 7: Update `OINInlineHandler.parse_inline_string()` in `src/oinsmiles/oin/inline.py`:
  - Remove `Chem.MolFromSmiles()` and `Chem.MolToSmiles()` calls.
  - Replace with: `smiles = SLOT_REGEX.sub('', oin_str); smiles = METAL_REGEX.sub('[M]', smiles); return smiles.strip('.')` (strip leading/trailing `.` separators).
  - Update docstring to note: "Returns non-canonical SMILES. @/@@ markers are preserved."
- [ ] Task 8: Wire into `XYZToSMILES.convert()` in `src/oinsmiles/core/translator.py`:
  - After `get_tmc_mol()`, call `Chem.SanitizeMol(mol)` then `CIPAssigner().assign_all(mol)`.
  - After `OINSanitizer.generate_robust_smiles()`, call `ChiralityRecoveryUtility().recover(mol)`.
  - Before serialization, call `PseudoAtomStrategy.strip_pseudo_atoms(mol)`.
- [ ] Task 9: Export `CIPAssigner`, `ChiralityRecoveryUtility`, `PseudoAtomStrategy` from `src/oinsmiles/__init__.py` if needed by downstream callers.

---

## 4. The Negative Space (Constraints)

* **DO NOT** call `CIPAssigner.assign_all()` after fragmentation — the metal context is required for correct CIP assignment of Zone A P/N atoms.
* **DO NOT** use Oganesson (Z=118) or `[Zz]` as the `PseudoAtomStrategy` sentinel — use `PSEUDO_ATOMIC_NUM = 0` (RDKit wildcard `*`).
* **DO NOT** allow `PseudoAtomStrategy` wildcard atoms (`*`) to appear in the final serialized OIN string — `strip_pseudo_atoms()` is mandatory before output.
* **DO NOT** modify `OINDiscreteAligner` — slot assignment logic is unchanged by this MiniPRD.
* **DO NOT** assert `@`/`@@` directly on Zone A atoms if the metal provides the 4th substituent — CIP assignment is only valid when the full coordination sphere is present (pre-fragmentation).
* **DO NOT** call `Chem.SanitizeMol()` inside `CIPAssigner.assign_all()` silently — if sanitization fails, the exception must propagate to the caller.

---

## 5. Integration Tests & Verification

* **Test 1 (Deterministic):** `CIPAssigner().assign_all(Chem.MolFromSmiles("N[C@@H](Cl)F"))` → returns mol; P atom at index... (use a known P-chiral SMILES for testing); `_CIPCode` property set on that atom.
* **Test 2 (Deterministic):** `OINInlineHandler().parse_inline_string("[Pt_SPL].[N@@H3]{0}.[Cl]{1}.[Cl]{2}.[N@@H3]{3}")` → returns SMILES containing `[N@@H3]` (not `[NH3]` — `@@` preserved).
* **Test 3 (Deterministic):** `PseudoAtomStrategy` fallback triggered when no `_CIPCode`: `strip_pseudo_atoms(mol)` returns mol with no atoms of `GetAtomicNum() == 0`.
* **Test 4 (Deterministic):** `CIPAssigner().assign_all(mol)` where mol is NOT sanitized → `Chem.SanitizeMol` exception propagates (does not silently no-op).
* **Test 5 (Novel):** P-chiral fixture XYZ → `XYZToSMILES.convert()` → OIN SMILES contains `@`/`@@` → [Candidate Artifact routing triggered — SMILES saved to `tests/candidate_outputs/pchiral_encoded.smi` for human review and RDKit CIP oracle verification before promotion to `tests/fixtures/`]
