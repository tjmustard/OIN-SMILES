# MiniPRD: Zone-A P Stereocenter Encoding — Encode Side (MiniPRD-A)
**Hypergraph Node ID:** atom_cip_assigner, atom_chirality_recovery
**Parent Node:** mod_core (`src/oinsmiles/core/chirality.py`)
**Parent SuperPRD:** `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` (v1.0.0)
**Execution order:** FIRST — establishes the `_OIN_CIPCode_LP` property contract MiniPRD-B consumes.

## 1. The Confidence Mandate
**Agent Instruction:** Before generating any plans or writing code, analyze this
document and output a Confidence Score (1-10). If the score is below 9, list
strictly the clarifying questions needed to reach 10.

Context an executor must load: `src/oinsmiles/core/chirality.py` (whole file),
`src/oinsmiles/utils/perception_tmc.py:30` (`TRANSITION_METALS_NUM`) and the
`recover()` call site at `perception_tmc.py:947-957`, SuperPRD §5.1 (dummy-copy recipe
+ data flow) and §6 (constraints).

## 2. Atomic User Stories
* **US-A1:** As the encode pipeline, I store a fragment-local (lone-pair) CIP
  label `_OIN_CIPCode_LP` on each P bonded to **exactly one** metal, computed by
  `rdCIPLabeler` on a dummy-metal copy while 3D is present, so the label
  survives fragmentation. Symmetric P → no property; ≥2-metal P → no property +
  `OINStereoWarning`; dummy-copy failure → no property + warning, `convert()`
  completes.
* **US-A2:** As the encode pipeline, `recover()` keeps and verify-and-flips a
  Zone-A P tag keyed on `_OIN_CIPCode_LP` (checked BEFORE any degree branch)
  instead of clearing it; multi-P fragments recompute CIP after each flip.
  DIPAMP emits `@`/`@@` on both `P{0}` and `P{1}`.
* **US-A3:** As a maintainer, `PseudoAtomStrategy` and all references (code,
  docstrings, comments, `architecture.yml` node + dangling edge) are deleted;
  the ≥4-neighbours-no-CIP clearing behaviour survives under a neutral name.
* **US-A4:** As the user, I get an `OINStereoWarning` (atom index in the message
  string) only on genuine same-convention conflict: `rdCIPLabeler` re-run **on
  the dummy-metal copy** vs the stored `_OIN_CIPCode_LP`. Never cross-convention.
  Skippable via a flag.
* **US-A5:** As the test suite, the DIPAMP OIN string is written to
  `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` as a Candidate Artifact —
  NOT promoted to a trusted golden without HITL sign-off (SuperPRD §9).

## 3. Implementation Plan (Task List)
- [ ] Task 1: Define `class OINStereoWarning(UserWarning)` in
      `core/chirality.py`; export it from the package namespace.
- [ ] Task 2: Import `TRANSITION_METALS_NUM` from `utils/perception_tmc.py` into
      `core/chirality.py` (import only — do NOT copy the list). If a circular
      import arises, move the constant to a new `core/constants.py` and import
      it from BOTH files (single source retained).
- [ ] Task 3: Write helper `_eligible_zone_a_p(mol) -> list[int]`: P atoms
      bonded to exactly one metal atom (predicate: atomic num in
      `TRANSITION_METALS_NUM`). P with ≥2 metal bonds → excluded +
      `warnings.warn(OINStereoWarning)` naming the atom index.
- [ ] Task 4: Write helper `_build_dummy_metal_copy(mol, p_idx) -> Chem.Mol | None`
      implementing the normative recipe: RWMol copy; metal atom → atomic num 0,
      formal charge 0, isotope 0; remove all metal–ligand bonds EXCEPT
      metal–P(p_idx); remove all bonds of any OTHER metal atoms;
      `Chem.SanitizeMol` — entire body in try/except; on any exception return
      `None` + `warnings.warn(OINStereoWarning, atom idx)`.
- [ ] Task 5: Write helper `_lp_cip_label(dummy_mol, p_idx) -> str | None` using
      `rdCIPLabeler.AssignCIPLabels` (guarded try/except; return None on
      error/no label). This is the ONLY label source for `_OIN_CIPCode_LP`.
- [ ] Task 6: In `CIPAssigner.assign_all()`: first clear any pre-existing
      `_OIN_CIPCode_LP` on all atoms (idempotence); then for each eligible P,
      build dummy copy → compute label → `SetProp('_OIN_CIPCode_LP', label)`;
      store nothing when label is None.
- [ ] Task 7: Add the diagnostic cross-check in `assign_all()` (skippable via
      `diagnostics: bool = True` parameter): re-run `_lp_cip_label` on the SAME
      dummy copy; on mismatch with the stored prop or None-where-tag-set,
      `warnings.warn(OINStereoWarning, f"... atom {idx} ...")`. Print (do not
      compare) the metal-present `_OIN_CIPCode` for HITL visibility.
- [ ] Task 8: Rewrite `ChiralityRecoveryUtility.recover()` branch order:
      (1) atom is P AND has `_OIN_CIPCode_LP` → keep chiral tag, recompute
      fragment-local CIP with `rdCIPLabeler` on the fragment mol, flip tag on
      mismatch — regardless of degree; (2) fall through to existing
      degree-keyed branches unchanged (Zone-A P/N without prop, all N → clear).
- [ ] Task 9: Multi-P handling in `recover()`: after any flip, re-run
      `rdCIPLabeler` on the fragment before evaluating the next tagged P;
      bound at 2 full passes over the tagged-P set.
- [ ] Task 10: Delete `class PseudoAtomStrategy` and `PSEUDO_ATOMIC_NUM` (verify
      unreferenced first with a repo grep); keep the ≥4-no-CIP clearing branch
      with a neutral comment; scrub the docstring reference at
      `tests/unit/test_stereo_roundtrip_diagnostics.py:181` and the `recover()`
      final-else comment.
- [ ] Task 11: Remove `atom_pseudo_atom_strategy` node AND the
      `atom_chirality_recovery.edges.depends_on → atom_pseudo_atom_strategy`
      edge from `spec/compiled/architecture.yml`.
- [ ] Task 12: Delete the duplicate fixture `tests/integration/Rh-RR-DIPAMP-Cl2.xyz`;
      canonical is `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz`.
- [ ] Task 13: Run `XYZToSMILES().convert()` on the canonical DIPAMP fixture;
      write the OIN string to `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt`
      with a provenance line (fixture path + SHA-256). Print the per-atom
      rdCIPLabeler table (both conventions, clearly headed) + emit a mol block
      for the HITL reviewer. Name the (R,R)-DIPAMP literature/CCDC reference in
      the sign-off note.
- [ ] Task 14: Tests — see §5. Run the full suite:
      `uv run python -m unittest discover tests`.

## 4. The Negative Space (Constraints)
* **DO NOT** duplicate the metal list — import `TRANSITION_METALS_NUM`
  (TD-005 lesson; a second list multiplies the stale-`is_metal` failure surface).
* **DO NOT** compare labels across conventions: the warning oracle runs on the
  dummy-metal copy ONLY; metal-present `_OIN_CIPCode` is print-only.
* **DO NOT** use legacy `AssignStereochemistry` `_CIPCode` values in any
  `_OIN_CIPCode_LP` computation, recomputation, or comparison — `rdCIPLabeler`
  end-to-end. (Legacy may still run for tag perception.)
* **DO NOT** let any dummy-copy or rdCIPLabeler exception escape `assign_all()`
  — degrade to store-nothing + `OINStereoWarning`.
* **DO NOT** derive a tag from 3D perception of a *trivalent* P — geometry↔tag
  work always uses the 4-coordinate dummy copy. (Graph-based CIP recompute from
  an existing tag on trivalent P in `recover()` is fine and required.)
* **DO NOT** store `_OIN_CIPCode_LP` on N, on ≥2-metal P, or on any atom where
  the label came back None.
* **DO NOT** emit `@`/`@@` on any nitrogen; keep clearing Zone-A N.
* **DO NOT** change `recover()` behaviour for atoms without `_OIN_CIPCode_LP`.
* **DO NOT** touch the OIN grammar, version string, or `oin/inline.py` parse code.
* **DO NOT** promote the DIPAMP candidate output to `tests/fixtures/` — HITL only.
* **DO NOT** regress carbon `@/@@` (TASK-10 set must stay green).

## 5. Integration Tests & Verification
* **Test 1 (Deterministic — property assignment):** `assign_all()` on the DIPAMP
  mol → both P atoms have `_OIN_CIPCode_LP` ∈ {R,S}; calling `assign_all()`
  twice yields identical properties (idempotence).
* **Test 2 (Deterministic — emitted string):** `XYZToSMILES().convert()` on
  `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` → OIN string contains `@` or `@@` on
  both `P{0}` and `P{1}` tokens.
* **Test 3 (Deterministic — negative controls):** BDPP + BDNN conversions are
  byte-identical to current goldens AND the emitted SMILES contains no
  `[P@`/`[P@@` token (explicit tag-absence assertion). Carbon-chirality TASK-10
  round-trips green.
* **Test 4 (Deterministic — raw parity, RISK-1):** on the DIPAMP fragment mol
  BEFORE `recover()`, compare rdCIPLabeler trivalent fragment CIP vs the stored
  `_OIN_CIPCode_LP` — assert equal, or mark the documented divergence (test
  makes divergence visible, not silent).
* **Test 5 (Deterministic — degradation, RISK-7):** CpM(PR₃)-type fixture
  (synthetic is fine): `convert()` completes, no `_OIN_CIPCode_LP` on the P if
  the dummy copy fails, output byte-identical to pre-feature output, exactly one
  `OINStereoWarning` naming the atom (use `warnings.catch_warnings(record=True)`).
* **Test 6 (Deterministic — bridging guard, B7):** synthetic mol with P bonded
  to two metals → no property + warning; conversion completes.
* **Test 7 (Deterministic — warning gate):** clean fixtures (DIPAMP, BDPP, BDNN,
  baseline complexes) pass under `-W error::OINStereoWarning`.
* **Test 8 (Deterministic — parse adjacency, B10):** `parse_inline_string` on a
  string containing `[P@]{0}` and `[P@@]{1>}` → slots, winding, and SMILES all
  correct (`_count_smiles_atoms_before` handles the bracket token).
* **Test 9 (Novel — Candidate Artifact):** `Rh-RR-DIPAMP-Cl2_oin.txt` written →
  **Candidate Artifact routing protocol triggered** — HITL review per SuperPRD
  §9 (lone-pair vs metal-present convention statement, rdCIPLabeler table,
  mol-block depiction, sign-off recorded in `spec/worklog/`).
