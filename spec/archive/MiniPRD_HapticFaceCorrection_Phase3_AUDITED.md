# MiniPRD: Haptic Face Correction (Stereo Phase 3)
**Hypergraph Node ID:** generation.molassembler_adapter.haptic_face_correction
**Parent Node:** SuperPRD_StereoPhase3_HapticFace (Stereo Phase 3 — Haptic Face Correction)
**Depends on:** MiniPRD_SignedCirculationHelper_Phase3 (must land first)

## 1. The Confidence Mandate
**Confidence Score**: 9/10.

The geometry (proper 180° in-plane C₂ rotation inverts `sign((v_star × v_next)·axis)`) is proven
sound by the Red Team. The residual 1 point is the empirical stability of the re-encoder's
heading-atom choice on the substituted fixture (R2), de-risked by the per-ring regex and the
explicit `xfail` fallback.

## 2. Atomic User Stories
* **US-001:** Input OIN eta-ring winding controls the generated haptic face; flipped `{n>}`↔`{n<}` yields a distinguishable, faithful 3D structure equal to the pinned golden.
* **US-003:** Ferrocene symmetric case handled honestly (documented skip; runtime symmetric ring is an identity no-op).
* **US-004:** Bridged ansa-metallocenes get winding-consistent placement where geometrically possible; conflicts are logged, not silently mis-placed.
* **US-005:** No behavior change on non-eta paths (byte-identical output; `test_winding_inertness.py` stays green).
* **US-006:** Correction operator proven chirality-preserving (`det +1` + CIP invariance).

## 3. Implementation Plan (Task List)
- [ ] Confirm `test_winding_inertness.py` is **green on `main`** before starting (regression attributability).
- [ ] In `_template_generate` (`:1034`): within each eta slot-group, read the single non-`None` `OINVector.winding` as the target marker; thread it as a scalar `winding` into `_stitch_eta_fragment`.
- [ ] Add an in-scope assertion that every in-scope eta geometry emits `OINVector`s carrying winding (guard the `oin_parser.py:495` template-gating hole); fail loudly if an in-scope shape loses it.
- [ ] **Multi-marker same-slot** (`c{0>}…c{0<}`): raise `ValueError` (canonical-form violation) — do not pick a winner.
- [ ] In `_stitch_eta_fragment`: after placement, compute measured circulation via `signed_circulation` (SMILES-ordered binding coords, `star_local_idx` tracking any analytic-path reorder at `:723-732`, metal→centroid axis).
- [ ] If measured char disagrees with target `winding`, build the correction `R` = proper 180° rotation about an in-plane axis through the centroid (⊥ metal→centroid). Axis = centroid→binding-atom[0] projected into the plane; **ε-guard**: if `‖projection‖ < 1e-6`, use binding-atom[1] then Gram–Schmidt against the axis.
- [ ] Assert `abs(np.linalg.det(R) - 1.0) < 1e-6` before applying (catches improper/reflection regressions at the source).
- [ ] Apply `R` to **all** fragment atoms; return a structured decision signal per ring (`fired` / `skipped` / `conflict` / `no-op`) instead of only `print(..., file=sys.stderr)`.
- [ ] **Zero-marker eta ring**: skip correction, record `no-op` in the signal (arbitrary-but-stable face).
- [ ] **Symmetric-ring detection**: identity no-op (no wasted rotation) at runtime.
- [ ] In `_stitch_multi_eta_fragment`: apply a coherent **whole-fragment** correction only when both rings disagree in the same sense; if they conflict (one wants a flip, the other does not), leave placement unchanged and record `conflict`.
- [ ] Pin `rdkit`/`numpy`/`scipy` versions used to bless the golden (§5.5 / `pyproject.toml`).

### Test tasks
- [ ] Create `tests/candidate_outputs/Ferrocene-halide-face_oin.txt` **pinned to the hand-verified string** `[Fe_LIN].Oc{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I.Oc{1}1[cH]{1}c{1}(I)c{1<}(Br)c{1}1Cl`; test fails if generation disagrees.
- [ ] Acceptance test on `Ferrocene-halide-face.xyz`: per-ring regex round-trip (anchored to a ring, not a specific heading atom); flip one ring's input → only that ring's char inverts.
- [ ] **Two-branch coverage**: assert one case where correction is skipped (pre-correction already matches) and one where it fires, via the instrumented decision signal — not merely by flipping the marker.
- [ ] Add a **chirality-witness fixture** (one eta ring, one pendant stereocenter); assert its CIP code is invariant across the correction.
- [ ] **Idempotency** test: double-applied correction (or generate→encode→generate) converges.
- [ ] Add a committed **bridged ansa-metallocene "before" baseline XYZ**; US-004 conflict path asserts no regression against it.
- [ ] Demote `test_haptic_face_winding` (ferrocene) to `@unittest.skip` with the symmetry-impossibility reason.
- [ ] R2 fallback (differ + geometric halide-sequence-reversal check) implemented as an explicit `xfail`/documented-skip with a logged reason — **never** a silent auto-downgrade.

## 4. The Negative Space (Constraints)
* **DO NOT** use a reflection (det −1) operator or negate a coordinate to "flip the face."
* **DO NOT** route the correction through `winding_by_slot` — use `OINVector.winding`.
* **DO NOT** duplicate the winding-sign math — call `signed_circulation`.
* **DO NOT** apply an independent per-ring correction inside a bridged multi-eta fragment.
* **DO NOT** change behavior on any non-eta path; output must be byte-identical.
* **DO NOT** silently auto-downgrade the exact assertion to the geometric fallback (a reflection also reverses the halide sequence, so it must never mask R1).
* **DO NOT** silently pick a winner for a multi-marker same-slot input — raise `ValueError`.
* **DO NOT** machine-emit-then-bless the golden — pin it to the reasoned §5.1 string.
* **DO NOT** overwrite `spec/active/Draft_PRD.md` or `RedTeam_Report_ZoneA_P_Encoding.md` (parallel Phase 4 session).

## 5. Integration Tests & Verification
* **Test 1 (Deterministic — golden):** generate from `Ferrocene-halide-face.xyz`'s OIN → re-encode → equals the pinned `Ferrocene-halide-face_oin.txt` string. Expected: exact match.
* **Test 2 (Deterministic — flip):** flip ring-0 input winding → ring-0 output char inverts, ring-1 unchanged. Expected: one char flips.
* **Test 3 (Deterministic — two branch):** decision signal shows one `skipped` and one `fired` across the two rings/cases. Expected: both branches exercised.
* **Test 4 (Deterministic — det +1):** `abs(det(R) − 1.0) < 1e-6`. Expected: proper rotation.
* **Test 5 (Deterministic — CIP invariance):** chirality-witness ring's CIP code unchanged across correction. Expected: invariant.
* **Test 6 (Deterministic — idempotency):** double correction converges. Expected: stable.
* **Test 7 (Deterministic — inertness):** non-eta XYZ output byte-unchanged; `test_winding_inertness.py` green. Expected: no diff.
* **Test 8 (Novel — Candidate Artifact):** `Ferrocene-halide-face.xyz` is a Candidate Artifact; if geometry quirks (0.95 Å O–H, Kabsch warning) destabilize encoding, regularize geometry and re-verify the pinned golden before blessing (routing protocol triggered).
