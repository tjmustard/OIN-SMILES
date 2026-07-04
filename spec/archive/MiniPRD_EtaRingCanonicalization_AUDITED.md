# MiniPRD: Eta-Ring Canonicalization (Fragment Order + Heading Atom)
**Hypergraph Node ID:** atom_oin_aligner
**Parent Node:** module_xyz_to_oin (XYZ→OIN encoder pipeline)
**Parent SuperPRD:** `spec/compiled/SuperPRD_EtaRingCanonicalization.md`
**Target file:** `src/oinsmiles/utils/oin_aligner.py` (`OINDiscreteAligner`)
**Executor tier:** Sonnet

## 1. The Confidence Mandate
**Agent Instruction:** Before writing code, output a Confidence Score (1-10). If below 9, list
the clarifying questions needed to reach 10.

**Pre-loaded confidence: 9/10.** Root causes and both blocking dispositions are pinned. The one
open empirical fact the executor MUST verify first (Task 1) is whether the lowest-canonical-rank
atom of each halide ring equals the atom the golden marks (`Oc{0<}` → O-carbon). If it does not,
STOP and report — do NOT re-pin the golden.

## 2. Atomic User Stories
* **US-001:** As the encoder, I want substituted-eta-ring fragments ordered by a
  heading-independent content-canonical key (scoped eta-only rank swap) so fragment order is
  invariant to xyz2mol arrival order, WITHOUT reordering any non-eta fragment.
* **US-002:** As the encoder, I want the heading atom of a substituted eta ring chosen as the
  lowest `Chem.CanonicalRankAtoms(breakTies=True)` atom (topology, not 3D orientation), with the
  `SYMMETRIC_LIGANDS` path unchanged and first-wins.
* **US-003:** As a maintainer, I want `test_haptic_face_golden_match` to flip from
  `expectedFailure` to a real byte-for-byte pass with `discover tests/unit` → skipped=3,
  expected failures=0, and every existing golden byte-identical.
* **US-004:** As a safety reviewer, I want proof the change is a pure winding-preserving
  relabeling: character start-invariant on the fixture ring, and a reflected ring still flips the
  marker; the R2 skip stays meaningful.

## 3. Implementation Plan (Task List)

**Phase 0 — Verify preconditions (do this BEFORE editing src):**
- [ ] Task 1: Write a throwaway probe: for each of the two halide rings in
  `tests/fixtures/Ferrocene-halide-face.xyz`, build the ring mol and compute
  `Chem.CanonicalRankAtoms(mol, breakTies=True)`; confirm the lowest-rank ring atom is the one
  the golden marks (`Oc{0<}` and `c{1<}(Br)`). If NOT equal → STOP and report; the heading rule
  or the golden is inconsistent — do not proceed and do not re-pin the golden.
- [ ] Task 2: Confirm `constituent_indices` (sorted local_idx of the ring atoms) tracks ring
  cyclic order for the fixture rings (so `signed_circulation`'s single-edge sign is
  start-invariant). Capture the evidence in the execution log.

**Phase 1 — RC2 heading (substituted eta only):**
- [ ] Task 3: In `_permute_and_serialize` (`:508-563`), before the geometric `best_idx` loop,
  add a content-canonical branch for eta groups NOT in `SYMMETRIC_LIGANDS`: obtain a fragment
  mol whose atom indexing maps to `constituent_indices` (RT-2: thread/rebuild with a
  round-trip-free index map; do NOT trust bare `MolFromSmiles(lig["smiles"])` indexing without
  a verified map), compute `Chem.CanonicalRankAtoms(mol, breakTies=True)`, and set
  `best_idx = the constituent atom with the lowest canonical rank`.
- [ ] Task 4: Fail-safe: if the mol/rank cannot be computed or mapped, fall back to the existing
  geometric `best_idx` for that fragment (RT-3). Keep the `SYMMETRIC_LIGANDS` override
  (`:565-597`) first-wins and untouched.

**Phase 2 — RC1 scoped eta-only rank swap:**
- [ ] Task 5: Add a heading-independent ring signature helper: canonical ring SMILES
  (`Chem.MolToSmiles`, order-invariant) computed via the SAME perception path that produces
  `lig["smiles"]` (RT-3), returning `None` on failure.
- [ ] Task 6: After slot/heading assignment, identify the set of **same-mass eta fragments** and
  the multiset of rank slots they occupy. Sort those fragments by
  `(canonical_ring_smiles, winding_sense_tiebreak, lowest_constituent_global_idx)` and reassign
  them to those same rank slots in sorted order. Winding sense enters the key ONLY to break ties
  between content-identical rings; `lowest_constituent_global_idx` is the final deterministic
  tiebreak (RT-5). Non-eta fragments and metal (rank 0) are never touched.
- [ ] Task 7: Retire the dead `base_sort_key` (`:250`, `"key"` at `:308`) — remove it or leave a
  one-line comment noting it is superseded by the scoped rank swap. (RT-1.)

**Phase 3 — Tests & acceptance:**
- [ ] Task 8: Remove `@unittest.expectedFailure` from `test_haptic_face_golden_match`
  (`tests/unit/test_stereo_roundtrip_diagnostics.py:469`); it must now pass byte-for-byte.
- [ ] Task 9: Add the four RT-5 hard tests (see §5).
- [ ] Task 10: Run acceptance commands (§5). Update `spec/worklog/NOTES.md` (Log entry +
  "HAPTIC CANON" task row → DONE). Do NOT git commit (leave to a separate step).

## 4. The Negative Space (Constraints)
* **DO NOT** modify `_determine_winding` or `oin.winding.signed_circulation` — the winding
  CHARACTER stays geometry-derived (FROZEN).
* **DO NOT** re-rank any non-eta fragment, change cross-species interleaving, or move the metal
  off rank 0. The RC1 swap permutes only same-mass eta fragments among the rank slots they
  already occupy.
* **DO NOT** touch or bypass the `SYMMETRIC_LIGANDS` forced-heading path (it stays first-wins).
* **DO NOT** let winding sense enter the RC1 key except as the content-identical tiebreak.
* **DO NOT** delete, weaken, un-skip, or auto-substitute
  `test_haptic_face_r2_geometric_fallback_never_auto_substituted`.
* **DO NOT** re-pin `tests/candidate_outputs/Ferrocene-halide-face_oin.txt` or any other golden.

## 5. Integration Tests & Verification

**Acceptance commands:**
- `uv run python -m unittest tests.unit.test_stereo_roundtrip_diagnostics -v`
- `uv run python -m unittest discover tests/unit` → **OK, skipped=3, expected failures=0**
- `uv run python -m unittest discover tests` → **OK**

* **Test 1 (Deterministic — golden, US-003):** `test_haptic_face_golden_match` un-`expectedFailure`d.
  Input: `Ferrocene-halide-face.xyz` → OIN(1) → generate → re-encode. Expected: byte-for-byte ==
  `tests/candidate_outputs/Ferrocene-halide-face_oin.txt`.
* **Test 2 (Non-eta inertness, RT-5):** encode cisplatin, transplatin, cis-PtCl2(en),
  fac/mer-Ir(ppy)3, PdCl2-RR-BDPP, PdCl2-RR-BDNN, PdCl2-R-BINAP → each byte-identical to its
  existing golden (guards R1: canonical `N` must NOT sort before `[Cl]`).
* **Test 3 (Symmetric-eta inertness, RT-5):** plain ferrocene and an ansa-metallocene golden →
  byte-identical (the `SYMMETRIC_LIGANDS` path must still win).
* **Test 4 (RC1 scoped-swap, RT-5):** an eta ring at a non-1 rank keeps every non-eta rank fixed
  (assert the non-eta slice of the output is byte-identical to the pre-fix output).
* **Test 5 (RC2 start-invariance + reflection, US-004/RT-4):** on a fixture halide ring, the
  winding character is identical for EVERY choice of star (`for star in ring_atoms: assert
  char(star) == char(ring_atoms[0])`); and a synthetically reflected ring still yields the
  flipped character (`assert char(reflected) != char(original)`). This is a LIVE assertion in
  addition to the kept R2 skip.
