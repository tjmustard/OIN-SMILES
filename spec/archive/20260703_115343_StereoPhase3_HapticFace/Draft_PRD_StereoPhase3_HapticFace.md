# Draft PRD: Stereo Phase 3 — Haptic Face Control (OIN→XYZ generation)

## Metadata
- **Project Name**: OIN-SMILES — Stereo Phase 3 (Haptic Face Control)
- **Version**: 0.1.0 (Draft — pre `/hyper-redteam`)
- **Status**: Draft (pre-RedTeam)
- **Owner**: Architect Agent (session 2026-07-03)
- **Roadmap**: `spec/worklog/ROADMAP-stereo.md` Phase 3
- **Depends on**: Stereo Phase 1 (winding plumbing) — DONE, `spec/archive/MiniPRD_WindingPlumbing_Phase1_AUDITED.md`
- **Note**: `spec/active/Draft_PRD.md` is occupied by the parallel Phase 4 (Zone-A P
  stereocenter) session; this Phase-3 draft uses a distinct filename to avoid clobbering it.
  Run `/hyper-redteam` against **this file** for Phase 3.

## 1. Introduction & Goals

### 1.1 Problem Statement
The OIN→XYZ generation direction silently ignores eta-ring winding. Phase 1 threaded
the V3.6 winding marker (`{n>}`/`{n<}`) as far as `ParsedOIN.winding_by_slot` and
`OINVector.winding`, but the 3D builder never consumes it: `_stitch_eta_fragment`
(`src/oinsmiles/generation/molassembler_adapter.py:584`) and
`_stitch_multi_eta_fragment` (`:151`) place eta rings by SVD centroid-plane alignment,
which fixes the ring normal toward/away from the metal but leaves the **circulation
sense of the ring atoms about the metal→centroid axis arbitrary** (whatever ETKDG's
embedding chirality produced). Consequently the prochiral face of a substituted Cp/arene
that coordinates the metal is not controllable from the input OIN. The Phase-0 diagnostic
(`tests/unit/test_stereo_roundtrip_diagnostics.py::test_haptic_face_winding`) captures this:
flipping a ring's input winding produces **byte-identical** 3D output.

### 1.2 Solution Overview
After eta placement, measure the signed circulation of the placed ring's binding atoms
(in fragment/SMILES order) about the metal→centroid axis; if its sign disagrees with the
winding threaded in from `OINVector.winding`, apply a **proper 180° rotation about an
in-plane axis through the ring centroid** (perpendicular to the metal→centroid axis) to the
whole fragment. This swaps which prochiral face coordinates the metal, reverses the measured
circulation, and — being a proper rotation (det +1) — preserves the chirality of any pendant
substituent. The signed-circulation sign computation is extracted into a single shared helper
so the encoder (`_determine_winding`) and the generation-side correction cannot drift.

### 1.3 Target Audience
Internal: `OIN3DGenerator.generate()` callers who supply an input OIN with eta rings and
expect the encoded winding to control the generated haptic face. No public-API surface change.

## 2. Confidence Mandate
**Confidence Score**: 9/10

All five interview branches were resolved with the user; the one residual uncertainty
(10 → 9) is the empirical stability of the re-encoder's heading-atom choice on the substituted
fixture, which the acceptance test's per-ring regex and documented weaker fallback already
de-risk.

**Clarifying Questions** (all resolved during interview):
- [x] Core operation is a face-swap that reverses circulation, superseding the roadmap's
  literal "mirror across the ring plane" wording. **Resolved:** accepted.
- [x] Winding source channel. **Resolved:** consume per-vector `OINVector.winding`, not
  `winding_by_slot` (the integer slot is not in hand at the placement site).
- [x] Multi-eta (bridged) scope. **Resolved:** full correction for single-eta; coherent
  whole-fragment reflection only for bridged multi-eta, with independent bridged-ring flip
  deferred to a future tether-dihedral-rotation mechanism.
- [x] Acceptance fixture. **Resolved:** ferrocene is symmetric → winding is not a geometric
  observable (unsatisfiable by geometry-only correction); re-target onto the user-provided
  `tests/integration/Ferrocene-halide-face.xyz` (each ring desymmetrized by H/OH/Cl/Br/I);
  demote the ferrocene test to a documented skip.
- [x] Correction operator. **Resolved:** proper 180° in-plane rotation (det +1, preserves
  pendant chirality), not a reflection.
- [x] Convention parity. **Resolved:** extract a shared `signed_circulation` helper as the
  single source of truth for both encode and generate sides.
- [x] Verification protocol. **Resolved:** Candidate-Artifact golden + faithful per-ring
  round-trip as primary assertion, with a documented weaker geometric fallback.

## 3. Scope

### 3.1 In-Scope
- Shared `signed_circulation(coords, star_local_idx, axis) -> '>'/'<'` helper (new module,
  e.g. `src/oinsmiles/oin/winding.py`); refactor `_determine_winding` to call it.
- Consume `OINVector.winding` at `_template_generate` and thread a scalar
  `winding: Optional[str]` into `_stitch_eta_fragment`; thread per-slot-group winding into
  `_stitch_multi_eta_fragment`.
- In `_stitch_eta_fragment`: after placement, if measured circulation disagrees with target
  winding, apply the proper 180° in-plane rotation to all fragment atoms.
- In `_stitch_multi_eta_fragment`: apply a coherent **whole-fragment** correction only when
  both rings disagree in the same sense; if they conflict, leave as-is and log the documented
  limitation.
- New acceptance test on `Ferrocene-halide-face.xyz` (faithful per-ring round-trip).
- Demote `test_haptic_face_winding` (ferrocene) to `@unittest.skip` with the
  symmetry-impossibility reason recorded.
- Candidate-Artifact golden: `tests/candidate_outputs/Ferrocene-halide-face_oin.txt`.

### 3.2 Out-of-Scope
- **Independent** face control of a single ring within a **bridged** multi-eta fragment
  (ansa-metallocene) — deferred to the future tether-dihedral-rotation mechanism (the 180°
  in-plane rotation is the free-ligand degenerate case of that same operation).
- Ligand P/N `@/@@` stereocenter encoding/enforcement (Phase 4 — parallel session).
- Any change to the metal-isomer / slot-ordering machinery (already correct — not a gap).
- Winding on non-eta or template-less (`NON`) geometries.
- Alternative structure-builder evaluation (Phase 4 decision).

## 4. User Stories (Atomic)

| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As a generation caller, I want the input OIN's eta-ring winding to control the generated haptic face, so a flipped `{n>}`↔`{n<}` yields a distinguishable, faithful 3D structure. | 1. On `Ferrocene-halide-face.xyz`, XYZ→OIN(1)→generate→re-encode reproduces each ring's input winding char (per-ring regex). 2. Flipping one ring's input winding inverts that ring's output char while the untouched ring is unchanged. | High |
| US-002 | As a maintainer, I want the encode and generate sides to share one winding-sign definition, so the correction cannot silently invert due to convention drift. | 1. `_determine_winding` and the placement correction both call `signed_circulation`. 2. A parity test asserts a placed geometry's helper-flip actually inverts the re-encoded marker. | High |
| US-003 | As a maintainer, I want the ferrocene symmetric case handled honestly, so the suite documents rather than fakes an unachievable pass. | 1. `test_haptic_face_winding` is `skip` with a recorded reason. 2. No non-physical winding metadata is injected to force a pass. | Medium |
| US-004 | As a generation caller of bridged ansa-metallocenes, I want winding-consistent placement where geometrically possible, so bridged output is not worse than today. | 1. Coherent whole-fragment correction applied when both rings agree. 2. Conflicting case logs the documented limitation and leaves placement unchanged (no regression vs current output). | Low |
| US-005 | As a maintainer, I want no behavior change on non-eta paths, so the change is inert outside haptic ligands. | 1. Full unit suite green. 2. Non-eta template/DG XYZ output byte-unchanged (pre/post diff harness). | High |

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Resolved Trade-offs

**Data flow.** `OIN3DGenerator.generate()` → `OINParser.parse()` → `ParsedOIN` (carries
`OINVector.winding` per binding vector + `winding_by_slot`) → `_template_generate` (`:1034`)
groups vectors per fragment, detects eta, and calls the stitch functions. The correction is
applied **inside** the stitch function, before positions are appended to `all_pos`.

**Winding channel (resolved).** Consume the per-vector `OINVector.winding`: within each eta
slot-group the single non-`None` value (heading atom) is the target marker. `winding_by_slot`
is keyed by integer slot, which is not available at the placement site (eta groups are keyed by
rounded vector-direction tuple), so it is not routed through the placement path.

**Correction operator (resolved).** Proper 180° rotation about an in-plane axis through the
ring centroid (⊥ metal→centroid axis). This reverses circulation about the axis *and* swaps the
coordinating face, det +1 → preserves pendant substituent chirality. A reflection (det −1) was
rejected: it would invert stereocenters inside substituents, colliding with Phase 4. The axis is
chosen deterministically (e.g. centroid→binding-atom[0] projected into the ring plane); its exact
choice only sets the final rotamer phase, which the existing post-placement ring-rotation
optimiser settles. Per user: this is the free-ligand degenerate case of a tether-dihedral flip —
same output.

**Optimiser ordering (derived, safe).** The post-placement optimiser (`:1170`) only *rotates*
each eta ring about the metal→centroid axis, which preserves circulation sign. A correction
applied earlier inside the stitch function therefore survives the optimiser untouched. Order:
correct-inside-stitch → append → optimiser spins harmlessly.

**Symmetry impossibility (resolved, drives fixture choice).** A regular unsubstituted Cp maps
onto itself as a point-set under the face-swap (its 5 vertices are a set-invariant); all ring
atoms are C and H is radial-symmetric, so the swapped structure is geometrically identical and
re-encodes to the identical marker. Hence winding is **not a geometric observable for symmetric
rings** — the roadmap's open question resolves to "V3.6 winding pins a physically meaningful
prochiral face **only for substituted rings**." Acceptance is therefore re-targeted to the
desymmetrized `Ferrocene-halide-face.xyz`, verified this session to encode as:
`[Fe_LIN].Oc{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I.Oc{1}1[cH]{1}c{1}(I)c{1<}(Br)c{1}1Cl`.

**Convention parity (resolved).** The sign computation `sign(cross(v_star, v_next) · slot_z)`
is extracted into `signed_circulation`, called by both `_determine_winding` (encode) and the
placement correction (generate), eliminating the sign-inversion bug class.

**Multi-eta bridged (resolved).** A bridged fragment cannot have one ring reflected without
tearing the Si bridge; only a whole-fragment operation is valid, which flips both rings
together. Applied only when both rings disagree in the same sense; conflicting windings log a
documented limitation. Independent bridged-ring control is deferred to a future
tether-dihedral-rotation mechanism (of which the 180° in-plane rotation is the free-ligand
degenerate case).

**Resolved Trade-offs Log:**
- **Issue:** Roadmap specifies "mirror across the ring plane." **Options:** reflect-across-ring-plane
  (a no-op for planar rings) vs reflection-through-axis-plane vs proper 180° in-plane rotation.
  **Resolution:** proper 180° in-plane rotation — only operator that reverses circulation, swaps
  the face, and preserves pendant chirality.
- **Issue:** Winding source `winding_by_slot` (per brief) vs `OINVector.winding`. **Resolution:**
  `OINVector.winding` — integer slot is unavailable at the placement site.
- **Issue:** Ferrocene acceptance test vs physical observability. **Resolution:** re-target to a
  desymmetrized fixture; demote ferrocene test to documented skip rather than fake a pass.
- **Issue:** Encode/generate sign convention drift. **Resolution:** single shared
  `signed_circulation` helper.

### 5.2 System Graph Blast Radius
Nodes affected (see `spec/compiled/architecture.yml`):
- `src/oinsmiles/generation/molassembler_adapter.py` — `_stitch_eta_fragment`,
  `_stitch_multi_eta_fragment`, `_template_generate` (thread + apply correction).
- `src/oinsmiles/utils/oin_aligner.py` — `_determine_winding` (refactor to call shared helper).
- **NEW** `src/oinsmiles/oin/winding.py` — `signed_circulation` shared helper.
- `tests/unit/test_stereo_roundtrip_diagnostics.py` — demote ferrocene test; add substituted test.
- **NEW** `tests/candidate_outputs/Ferrocene-halide-face_oin.txt` — Candidate-Artifact golden.
- `tests/integration/Ferrocene-halide-face.xyz` — acceptance fixture (already provided by user).

### 5.3 Execution Checklist (MiniPRDs)
To be generated by `/hyper-resolve` after RedTeam:
- [ ] `spec/compiled/MiniPRD_SignedCirculationHelper_Phase3.md` — shared helper + `_determine_winding` refactor + parity test.
- [ ] `spec/compiled/MiniPRD_HapticFaceCorrection_Phase3.md` — thread winding, apply 180° correction in stitch functions, multi-eta coherent case, acceptance test + ferrocene skip.

(Single combined MiniPRD acceptable if RedTeam finds the split unnecessary.)

### 5.4 API Contracts / Schema
```python
# NEW: src/oinsmiles/oin/winding.py
def signed_circulation(
    coords: np.ndarray,      # (n,3) placed binding-atom positions, fragment/SMILES order
    star_local_idx: int,     # index (into coords) of the heading/star atom
    axis: np.ndarray,        # metal→centroid outward axis (need not be unit)
) -> str:                    # '>' if cross(v_star, v_next)·axis >= 0 else '<'
    ...

# CHANGED signatures (internal only; sole caller is _template_generate)
def _stitch_eta_fragment(
    frag_smiles, binding_idxs, slot_unit, metal_sym,
    winding: Optional[str] = None,   # NEW, default None = no correction (back-compat)
) -> tuple[np.ndarray, list[str], "Chem.Mol | None"] | None: ...

def _stitch_multi_eta_fragment(
    frag_smiles, vectors, metal_sym,
    # per-slot winding read from vectors' .winding; no new positional param required
) -> tuple[np.ndarray, list[str], "Chem.Mol | None"] | None: ...
```

### 5.5 Dependencies
- Existing only: `numpy`, `scipy.spatial.transform.Rotation`, `rdkit`. No new third-party deps.

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** use a reflection (improper) operator — it inverts pendant substituent chirality.
- **DO NOT** interpret "mirror across the ring plane" literally — it is a no-op for planar rings.
- **DO NOT** inject winding as non-physical metadata into the XYZ to force the ferrocene test.
- **DO NOT** route the placement correction through `winding_by_slot` (integer slot unavailable
  at the site); use `OINVector.winding`.
- **DO NOT** duplicate the winding-sign math in the adapter — call `signed_circulation`.
- **DO NOT** apply an independent per-ring correction inside a bridged multi-eta fragment.
- **DO NOT** change behavior on any non-eta path; output there must be byte-identical.
- **DO NOT** overwrite `spec/active/Draft_PRD.md` (owned by the parallel Phase 4 session).

## 7. Risks & Mitigation
- **R1 — Sign inversion (correction flips the wrong way).** → Shared `signed_circulation` helper
  + parity test asserting a placed-geometry flip actually inverts the re-encoded marker.
- **R2 — Re-encoder heading-atom choice unstable on the substituted fixture** (breaks exact
  round-trip). → Per-ring regex assertion (robust to heading identity); documented weaker
  fallback (differ + geometric halide-sequence-reversal check).
- **R3 — Optimiser clobbers the correction.** → Derived safe (axis-rotation preserves circulation
  sign); guarded by the byte-diff harness and acceptance test.
- **R4 — Fixture geometry quirks** (0.95 Å O–H, poorly-defined Kabsch warning). → Fixture is a
  Candidate Artifact; if quirks destabilize encoding, regularize geometry and re-verify golden.
- **R5 — Multi-eta conflicting windings silently mis-placed.** → Explicit logged limitation; no
  worse than current output; covered by US-004.

## 8. Success Metrics
- New substituted-fixture acceptance test passes as a hard assert (faithful per-ring round-trip).
- `test_haptic_face_winding` (ferrocene) present as a documented skip.
- Parity test (US-002) passes.
- Full `uv run python -m unittest discover tests/unit` green; non-eta XYZ output byte-unchanged.
- `signed_circulation` is the sole winding-sign definition (grep: no duplicated cross/dot sign
  logic in `molassembler_adapter.py`).
