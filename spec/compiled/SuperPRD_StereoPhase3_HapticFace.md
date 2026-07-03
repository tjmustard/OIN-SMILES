# SuperPRD — Stereo Roadmap Phase 3: "Control the Face" (Haptic Face Correction)

## Metadata
- **Project Name**: OIN-SMILES — Generation-side stereochemistry, Phase 3
- **Version**: 1.0.0 (Compiled — post `/hyper-resolve`)
- **Status**: Ready for `/hyper-execute`
- **Owner**: Architect Agent / Thomas Mustard
- **Roadmap ref**: `spec/worklog/ROADMAP-stereo.md` § Phase 3
- **Depends on**: Stereo Phase 1 (winding plumbing) — DONE, `spec/archive/MiniPRD_WindingPlumbing_Phase1_AUDITED.md`
- **Provenance**: `Draft_PRD_StereoPhase3_HapticFace.md` (v0.1.0) → `RedTeam_Report_StereoPhase3_HapticFace.md` (2026-07-03) → this compilation
- **Child MiniPRDs**:
  - `spec/compiled/MiniPRD_SignedCirculationHelper_Phase3.md`
  - `spec/compiled/MiniPRD_HapticFaceCorrection_Phase3.md`
- **Sibling (parallel, DO NOT clobber)**: Phase 4 Zone-A P-stereocenter session
  (`spec/active/Draft_PRD.md`, `spec/active/RedTeam_Report_ZoneA_P_Encoding.md`)

---

## 1. Introduction & Goals

### 1.1 Problem Statement
The OIN→XYZ generation direction silently ignores eta-ring **winding**. Phase 1 threaded the
V3.6 winding marker (`{n>}`/`{n<}`) as far as `ParsedOIN.winding_by_slot` and
`OINVector.winding`, but the 3D builder never consumes it. `_stitch_eta_fragment`
(`molassembler_adapter.py:584`) and `_stitch_multi_eta_fragment` (`:151`) place eta rings by
SVD centroid-plane alignment, which fixes the ring normal toward/away from the metal but leaves
the **circulation sense of the ring atoms about the metal→centroid axis arbitrary** (whatever
ETKDG's seed-42 embedding produced). The prochiral face of a substituted Cp/arene that
coordinates the metal is therefore not controllable from the input OIN. The Phase-0 diagnostic
(`test_stereo_roundtrip_diagnostics.py::test_haptic_face_winding`) proves it: flipping a ring's
input winding produces **byte-identical** 3D output.

### 1.2 Solution Overview
After eta placement, measure the signed circulation of the placed ring's binding atoms (in
fragment/SMILES order) about the metal→centroid axis; if its sign disagrees with the winding
threaded in from `OINVector.winding`, apply a **proper 180° rotation about an in-plane axis
through the ring centroid** (⊥ the metal→centroid axis) to the whole fragment. This swaps which
prochiral face coordinates the metal, reverses the measured circulation, and — being a proper
rotation (det +1) — preserves the chirality of any pendant substituent. The signed-circulation
computation is extracted into a single shared helper (`signed_circulation`) so the encoder
(`_determine_winding`) and the generation-side correction cannot drift.

### 1.3 Target Audience
Internal: `OIN3DGenerator.generate()` callers who supply an input OIN with eta rings and expect
the encoded winding to control the generated haptic face. No public-API surface change.

---

## 2. Confidence Mandate
**Confidence Score**: 9/10 (post-resolution)

The Red Team correctly downgraded the draft's 9/10 to 7/10 pending three load-bearing items
(helper contract, `det +1` verifiability, golden origin). All three are now resolved with
full-closure decisions below, restoring 9/10. The residual 1 point is the empirical stability of
the re-encoder's heading-atom choice on the substituted fixture (R2), de-risked by the per-ring
regex assertion and the documented weaker fallback (now an explicit `xfail`, never a silent
downgrade).

---

## 3. Resolved Decisions (Red Team → Resolution)

All eight ranked Red Team findings carry a documented decision.

| # | Finding | Resolution |
| :-- | :-- | :-- |
| 1 | `signed_circulation` centroid/ordering/axis contract underspecified (highest risk). | **Full closure.** The helper (a) centers `coords` on their own centroid **internally**; (b) requires `coords` in exact SMILES/fragment order, and `_determine_winding` is refactored to feed that same order; (c) hosts the `n < 3 → '>'` default **inside** the helper so both sides inherit it; (d) fixes `axis` orientation to metal→centroid (outward), with a unit test on an antiparallel/inward-Z slot asserting parity. Backed by a **fixture-independent behavioral parity test** (synthetic labeled pentagon fed through both call paths → identical char; then 180° rotation → char flips). |
| 2 | Chosen fixture has no stereocenter → a reflecting (det −1) bug ships green, detonates in Phase 4. | **det assert + chirality witness.** Add a runtime `abs(det(R) − 1.0) < 1e-6` assertion on the correction transform. Add a minimal **chirality-witness fixture** (one eta ring bearing a single pendant stereocenter) whose CIP code must be invariant across the correction. In-scope for *this operator's* correctness even though Phase 4 encoding is not. |
| 3 | Fixed ETKDG seed may leave one correction branch uncovered. | **Instrumented two-branch coverage.** The correction reports whether it fired (metadata flag). The acceptance test asserts one case where pre-correction circulation already matches target (rotation **skipped**) and one where it does not (rotation **fired**) — verified by instrumentation, not merely by flipping the input marker. |
| 4 | Golden does not exist yet; risks self-blessing a wrong artifact. | **Pin to hand-verified string.** `Ferrocene-halide-face_oin.txt` is set to the §5.1 hand-reasoned string and the test **fails** if generation disagrees. The golden is a reasoned spec artifact, never machine-emitted-then-blessed. |
| 5 | US-004 "conflict" undefined; "no regression" baseline unpinned. | **Concrete definition + committed baseline.** "Conflict" ≝ one ring's measured circulation wants a flip while the other's does not. Commit a "before" baseline XYZ for a bridged ansa-metallocene so "no regression" is falsifiable. |
| 6 | No observability of correction decisions (stderr prints only). | **Structured signal.** Replace `print(…, file=sys.stderr)` with an inspectable return-metadata signal reporting `fired` / `skipped` / `conflict` / `no-op` per ring, so silent mis-placement is detectable. |
| 7 | Degenerate in-plane axis has no specified guard. | **ε-guard + deterministic fallback.** If `‖centroid→binding-atom[0] projected into plane‖ < 1e-6`, fall back to binding-atom[1], then Gram–Schmidt against the axis. |
| 8 | Winding channel gated on template branch (`oin_parser.py:495`). | **Assert emission.** Add an in-scope assertion that every in-scope eta geometry emits `OINVector`s carrying winding; if an in-scope shape ever loses it, fail loudly rather than silently no-op. |

**Additional edge-case / NFR decisions (Phase 2):**
- **Zero-marker eta ring** (legacy V3.0–V3.5 or hand-authored, no marker): documented no-op —
  correction skipped and logged via the structured signal; face is arbitrary-but-stable.
- **Multi-marker same slot** (`c{0>}…c{0<}`): **raise `ValueError`** (canonical-form violation);
  never silently pick whichever the iterator hits first.
- **Idempotency**: applying the correction twice (or round-tripping generate→encode→generate)
  must converge, not oscillate — asserted.
- **Symmetric-ring runtime no-op**: on a detected-symmetric ring the correction is **identity**
  (not a wasted rotation), matching the "winding is not a geometric observable for symmetric
  rings" resolution, so runtime symmetric input is inert.
- **Version pin (NFR)**: pin `rdkit` / `numpy` / `scipy` in §5.5 — the entire diagnostic rests on
  seed-42 embedding determinism.

---

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As a generation caller, I want the input OIN's eta-ring winding to control the generated haptic face, so a flipped `{n>}`↔`{n<}` yields a distinguishable, faithful 3D structure. | 1. On `Ferrocene-halide-face.xyz`, XYZ→OIN(1)→generate→re-encode reproduces each ring's input winding char (per-ring regex, anchored to a ring not a specific heading atom). 2. Flipping one ring's input winding inverts that ring's output char while the untouched ring is unchanged. 3. Generation output equals the §5.1 hand-verified golden string. | High |
| US-002 | As a maintainer, I want encode and generate to share one winding-sign definition, so the correction cannot silently invert due to convention drift. | 1. Both `_determine_winding` and the placement correction call `signed_circulation`. 2. **Behavioral parity** — a synthetic placed ring's helper char equals the encoder's char on the identical coordinate array (fixture-independent). 3. **Operator invert** — applying the 180° rotation flips the re-encoded char (instrumented). | High |
| US-003 | As a maintainer, I want the ferrocene symmetric case handled honestly, so the suite documents rather than fakes an unachievable pass. | 1. `test_haptic_face_winding` is `skip` with a recorded symmetry-impossibility reason. 2. No non-physical winding metadata is injected to force a pass. 3. A detected-symmetric ring at runtime is an identity no-op (no wasted rotation). | Medium |
| US-004 | As a generation caller of bridged ansa-metallocenes, I want winding-consistent placement where geometrically possible, so bridged output is not worse than today. | 1. Coherent whole-fragment correction applied when both rings agree. 2. "Conflict" (one ring wants a flip, the other does not) logs the documented limitation via the structured signal and leaves placement unchanged. 3. A committed "before" baseline XYZ pins "no regression" as falsifiable. | Low |
| US-005 | As a maintainer, I want no behavior change on non-eta paths, so the change is inert outside haptic ligands. | 1. Full unit suite green, **including the pre-existing `test_winding_inertness.py`** (confirmed green on `main` first, so any regression is attributable). 2. Non-eta template/DG XYZ output byte-unchanged (pre/post diff harness). | High |
| US-006 | As a maintainer, I want the correction operator proven chirality-preserving, so a reflection bug cannot reach Phase 4. | 1. `abs(det(R) − 1.0) < 1e-6` asserted on the correction transform. 2. A chirality-witness eta ring's CIP code is invariant across the correction. | High |

---

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Data Flow
`OIN3DGenerator.generate()` → `OINParser.parse()` → `ParsedOIN` (carries `OINVector.winding`
per binding vector + `winding_by_slot`) → `_template_generate` (`:1034`) groups vectors per
fragment, detects eta, and calls the stitch functions. The correction is applied **inside** the
stitch function, before positions are appended to `all_pos`. The post-placement optimiser
(`:1170`) rotates only about `slot_u` (`_Rot.from_rotvec(slot_u * angle)`), which preserves
circulation sign — so a correction applied earlier survives untouched.

**Winding channel (resolved).** Consume per-vector `OINVector.winding`: within each eta
slot-group the single non-`None` value (heading atom) is the target marker. `winding_by_slot` is
keyed by integer slot, unavailable at the placement site (eta groups keyed by rounded
vector-direction tuple), so it is not routed through the placement path.

**Correction operator (resolved).** Proper 180° rotation about an in-plane axis through the ring
centroid (⊥ metal→centroid axis). Reverses circulation about the axis, swaps the coordinating
face, det +1 → preserves pendant substituent chirality. A reflection (det −1) is rejected: it
would invert stereocenters inside substituents, colliding with Phase 4. The axis is chosen
deterministically (centroid→binding-atom[0] projected into the ring plane; ε-guarded fallback to
binding-atom[1] + Gram–Schmidt when the projection is ill-conditioned). Its exact choice only
sets the final rotamer phase, which the existing post-placement ring-rotation optimiser settles.

**Hand-verified golden (Candidate Artifact, pinned).**
`Ferrocene-halide-face.xyz` encodes as:
```
[Fe_LIN].Oc{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I.Oc{1}1[cH]{1}c{1}(I)c{1<}(Br)c{1}1Cl
```
`tests/candidate_outputs/Ferrocene-halide-face_oin.txt` is pinned to this exact string; the
acceptance test fails if generation disagrees.

### 5.2 System Graph Blast Radius (confirmed contained)
- **NEW** `src/oinsmiles/oin/winding.py` — `signed_circulation` shared helper (centers
  internally, hosts `n<3` default, fixed axis convention).
- `src/oinsmiles/utils/oin_aligner.py` — `_determine_winding` (sole caller `oin_aligner.py:619`)
  refactored to call the shared helper with SMILES-ordered coords.
- `src/oinsmiles/generation/molassembler_adapter.py` — `_stitch_eta_fragment` (sole caller
  `_template_generate` at `:1084`, back-compat via `winding=None` default), `_stitch_multi_eta_fragment`,
  `_template_generate` (thread + apply correction + structured signal).
- `tests/unit/test_stereo_roundtrip_diagnostics.py` — demote ferrocene test to skip; add
  substituted acceptance test, two-branch coverage, behavioral parity, `det +1`, CIP-invariance,
  idempotency.
- **NEW** `tests/candidate_outputs/Ferrocene-halide-face_oin.txt` — pinned golden.
- `tests/integration/Ferrocene-halide-face.xyz` — acceptance fixture (user-provided).
- **NEW** chirality-witness fixture (one eta ring, one pendant stereocenter).
- **NEW** bridged ansa-metallocene "before" baseline XYZ (US-004 no-regression pin).
- Pre-existing `tests/unit/test_winding_inertness.py` must stay green (US-005 inertness contract).

### 5.3 API Contracts / Schema
```python
# NEW: src/oinsmiles/oin/winding.py
def signed_circulation(
    coords: np.ndarray,      # (n,3) placed binding-atom positions, EXACT SMILES/fragment order
    star_local_idx: int,     # index (into coords) of the heading/star atom
    axis: np.ndarray,        # metal→centroid OUTWARD axis (need not be unit)
) -> str:
    """Single source of truth for winding sign, called by BOTH encode and generate.

    Contract (Red Team Finding #1, resolved):
      - coords are centered on their OWN centroid INSIDE this function.
      - coords MUST be in SMILES/fragment order; v_next is the cyclic +1 neighbour.
      - n < 3 returns '>' (default lives HERE so both call sites inherit it).
      - axis convention is metal→centroid outward; encode-side slot_z feed must match
        (asserted by an antiparallel-slot unit test).
    Returns '>' if cross(v_star, v_next)·axis >= 0 else '<'.
    """

# CHANGED signatures (internal only; sole caller is _template_generate)
def _stitch_eta_fragment(
    frag_smiles, binding_idxs, slot_unit, metal_sym,
    winding: Optional[str] = None,   # NEW, default None = no correction (back-compat)
) -> tuple[np.ndarray, list[str], "Chem.Mol | None"] | None: ...

def _stitch_multi_eta_fragment(
    frag_smiles, vectors, metal_sym,
    # per-slot winding read from vectors' .winding; whole-fragment coherent correction only
) -> tuple[np.ndarray, list[str], "Chem.Mol | None"] | None: ...
```

**Correction-decision signal (Red Team #6, R6).** The stitch functions report per-ring outcome
(`fired` / `skipped` / `conflict` / `no-op`) as inspectable return metadata rather than only
`print(..., file=sys.stderr)`, so US-004 conflicts and zero-marker no-ops are detectable.

### 5.4 Dependencies (NFR — version pin)
- Existing only: `numpy`, `scipy.spatial.transform.Rotation`, `rdkit`. **No new third-party deps.**
- **Pin** the `rdkit` / `numpy` / `scipy` versions used to bless the seed-42 golden in §5.5 of the
  MiniPRD / `pyproject.toml`; record them so a future embedding-handedness change is a loud test
  failure, not a silent wrong-face golden.

---

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** use a reflection (improper, det −1) operator — it inverts pendant substituent chirality.
- **DO NOT** compute `signed_circulation` on raw, un-centered coords or in a non-SMILES order —
  centering and ordering are part of the helper's contract; violating them recreates the drift bug.
- **DO NOT** put the `n<3` default in only one call site — it lives inside the helper.
- **DO NOT** interpret "mirror across the ring plane" literally — it is a no-op for planar rings.
- **DO NOT** inject non-physical winding metadata into XYZ to force the ferrocene test.
- **DO NOT** silently pick a winner for a multi-marker same-slot input — raise `ValueError`.
- **DO NOT** silently auto-downgrade the exact assertion to the geometric fallback — the fallback
  is an explicit `xfail`/documented-skip with a logged reason (a reflection would also satisfy a
  halide-sequence-reversal check, so it must never mask R1).
- **DO NOT** route the placement correction through `winding_by_slot` — use `OINVector.winding`.
- **DO NOT** duplicate the winding-sign math in the adapter — call `signed_circulation`.
- **DO NOT** apply an independent per-ring correction inside a bridged multi-eta fragment.
- **DO NOT** change behavior on any non-eta path; output there must be byte-identical.
- **DO NOT** overwrite `spec/active/Draft_PRD.md` or `RedTeam_Report_ZoneA_P_Encoding.md`
  (owned by the parallel Phase 4 session).

---

## 7. Out-of-Scope
- **Independent** face control of a single ring within a **bridged** multi-eta fragment
  (ansa-metallocene) — deferred to the future tether-dihedral-rotation mechanism (the 180°
  in-plane rotation is that operation's free-ligand degenerate case).
- Ligand P/N `@/@@` stereocenter encoding/enforcement (Phase 4 — parallel session). Note: the
  chirality-witness fixture here tests *this operator's* CIP invariance, not Phase 4 encoding.
- Any change to the metal-isomer / slot-ordering machinery (already correct — not a gap).
- Winding on non-eta or template-less (`NON`) geometries beyond the emission assertion.
- Alternative structure-builder evaluation (Phase 4 decision).

---

## 8. Risks & Mitigation
- **R1 — Sign inversion.** → Shared `signed_circulation` (centered, ordered, single `n<3`
  default) + behavioral parity test on real placed coords (not a synthetic ideal only).
- **R2 — Re-encoder heading-atom instability on the substituted fixture.** → Per-ring regex
  (robust to heading identity); documented weaker fallback as an explicit `xfail`, never a silent
  auto-downgrade.
- **R3 — Optimiser clobbers the correction.** → Derived safe (axis-rotation preserves circulation
  sign, confirmed in code at `:1170`); guarded by byte-diff harness + acceptance test.
- **R4 — Fixture geometry quirks** (0.95 Å O–H, poor Kabsch warning). → Fixture is a Candidate
  Artifact; golden pinned to the reasoned §5.1 string, re-verified deterministic before blessing.
- **R5 — Multi-eta conflicting windings mis-placed.** → Explicit logged limitation via structured
  signal; committed baseline pins no-regression; covered by US-004.
- **R6 — Observability.** → Structured, inspectable correction-decision signal replaces stderr
  prints, so silent no-op / conflict paths are visible.
- **R7 — Degenerate in-plane axis.** → ε-guard + deterministic fallback axis (binding-atom[1] +
  Gram–Schmidt).
- **R8 — Reflection ships undetected on a stereocenter-free fixture.** → `det(R) ≈ +1` runtime
  assertion + chirality-witness ring CIP invariance.

---

## 9. Success Metrics
- New substituted-fixture acceptance test passes as a hard assert (faithful per-ring round-trip),
  and generation equals the §5.1 hand-verified golden.
- `test_haptic_face_winding` (ferrocene) present as a documented skip.
- **Behavioral parity** (US-002.2): a placed ring's helper char equals the encoder's char on the
  identical coordinate array.
- **Operator invert** (US-002.3): the correction flips the re-encoded char (instrumented, both
  branches — skipped and fired).
- **`det(R) ≈ +1`** on the correction transform (US-006.1).
- **CIP invariance** on the chirality-witness ring (US-006.2).
- **Idempotency**: double-applied correction converges.
- Full `uv run python -m unittest discover tests/unit` green (incl. `test_winding_inertness.py`);
  non-eta XYZ output byte-unchanged.
- `signed_circulation` is the sole winding-sign definition (grep: no duplicated cross/dot sign
  logic in `molassembler_adapter.py`) — a **textual** check backstopped by the behavioral parity
  metric above.

---

## 10. Execution Plan (Child MiniPRDs)
- [ ] `spec/compiled/MiniPRD_SignedCirculationHelper_Phase3.md` — shared helper (centered/ordered
  contract, `n<3` default, fixed axis convention) + `_determine_winding` refactor + antiparallel-slot
  unit test + fixture-independent behavioral parity test.
- [ ] `spec/compiled/MiniPRD_HapticFaceCorrection_Phase3.md` — thread winding, apply 180° in-plane
  correction in stitch functions with `det +1` assertion, structured decision signal, degenerate-axis
  guard, multi-eta coherent/conflict handling, edge-case contract (zero-marker no-op, multi-marker
  ValueError, symmetric-ring identity, idempotency), acceptance test + two-branch coverage + pinned
  golden + chirality-witness + ferrocene skip.
