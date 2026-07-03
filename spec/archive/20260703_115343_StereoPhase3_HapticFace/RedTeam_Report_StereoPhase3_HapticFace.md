# Red Team Report: Stereo Phase 3 — Haptic Face Control

- **Target PRD**: `spec/active/Draft_PRD_StereoPhase3_HapticFace.md` (v0.1.0)
- **Reviewer**: Red Team Agent (session 2026-07-03)
- **Blast-radius source**: `spec/compiled/architecture.yml`, live source under `src/oinsmiles/`
- **Verdict**: The design is unusually well-resolved for a pre-RedTeam draft (the geometry
  of the 180° in-plane correction is sound — a C₂ rotation about an in-plane axis ⊥ the
  metal→centroid axis does provably invert `sign((v_star × v_next)·axis)`). The confidence
  claim of **9/10 is too high**: at least three load-bearing details are underspecified in a
  way that reintroduces the exact sign-drift bug class the PRD claims to eliminate, and the
  chosen acceptance fixture cannot exercise the one property (`det +1` preserves pendant
  chirality) that justifies the correction operator over a reflection.

> Grounding notes carried into the analysis below (verified against source this session):
> - `_stitch_eta_fragment` sole caller is `_template_generate` (`molassembler_adapter.py:1084`) — back-compat claim holds.
> - `_determine_winding` sole caller is `oin_aligner.py:619` — refactor blast radius is contained.
> - The `:1170` optimiser rotates only about `slot_u` (`_Rot.from_rotvec(slot_u * angle)`) — "optimiser preserves circulation" (R3) is confirmed by code.
> - `OINVector` objects are created **only inside the template branch** (`oin_parser.py:495`, gated on `tmpl_vectors is not None and sa.slot < len(tmpl_vectors)`).
> - `tests/candidate_outputs/Ferrocene-halide-face_oin.txt` **does not exist yet**.
> - `tests/unit/test_winding_inertness.py` already exists (relevant to US-005).

---

## Section 1 (Introduction & Goals) Analysis

* **Clarifying Questions:**
  - The problem statement says ETKDG's "embedding chirality" makes the pre-correction
    circulation arbitrary, yet `_stitch_eta_fragment` pins `randomSeed = 42`. With a fixed
    seed the pre-correction sign is **deterministic per fixture**, not random. Does that mean
    for the acceptance fixture the correction branch might *never fire* for one of the two
    windings (ETKDG happens to already produce it), so the flip test only exercises the
    rotation in one direction? What guarantees both the "correction-applied" and
    "correction-skipped" branches are covered?
  - The analytic-fallback path (`_analytic_ring_geometry`) builds a ring with a **fixed
    rotational phase** (`positions.append(...cos/sin...)`). Its circulation sign is therefore
    a constant of the code, not of the input. Is the correction expected to fire on the
    analytic path too, and has the constant sign been measured (not assumed)?
  - "Byte-identical 3D output" is cited as the diagnostic. Byte-identity is sensitive to
    float formatting; is the diagnostic asserting geometric equality or literal string
    equality, and will the new test use the same notion?

* **What-If Scenarios:**
  - **Seed masks a branch.** Seed 42 yields the `>` face for both rings pre-correction. The
    `{n<}` input then always triggers the rotation and `{n>}` never does. A sign bug in the
    *skip* path (or an off-by-one in "which direction is already correct") ships green because
    the test only ever compares against the one branch ETKDG produced.
  - **Embedding instability across RDKit versions.** A future RDKit bump changes the seed-42
    embedding handedness. The golden `Ferrocene-halide-face_oin.txt` was blessed under the old
    handedness; the correction now fires on the opposite input, and the golden silently
    encodes the *wrong* face while still "round-tripping."

* **Points for Improvement:**
  - Add an explicit **two-branch coverage requirement**: one assertion where the pre-correction
    circulation already matches the target (no rotation) and one where it does not (rotation
    fires), verified by instrumenting whether the rotation was applied — not merely by flipping
    the input marker.
  - State the RDKit/numpy/scipy version pin (currently absent from §5.5) as an NFR, since the
    entire diagnostic rests on embedding determinism.

---

## Section 2 (Confidence Mandate) Analysis

* **Clarifying Questions:**
  - The single stated residual (heading-atom stability, R2) is real, but why is the
    **frame/centroid contract of `signed_circulation`** (see §5 below) not counted as a
    residual? The encoder subtracts the centroid before the cross product
    (`oin_aligner.py:676-677`), while the §5.4 API takes raw `coords` + `axis` and no centroid.
    If the shared helper does not center internally, encode and generate compute *different*
    cross products and the "single source of truth" still drifts.
  - What is the evidence basis for 9/10 vs 7/10? No prototype of the correction is claimed to
    have been run end-to-end on the fixture; the encoded golden string in §5.1 is asserted
    "verified this session," but generation-side correction is not.

* **What-If Scenarios:**
  - The confidence number is used by `/hyper-resolve` to decide MiniPRD granularity. A 9/10
    that is actually 7/10 leads resolve to under-scope tests and skip the parity/idempotency
    guards, shipping the sign-drift class it claims to have closed.

* **Points for Improvement:**
  - Downgrade to **7/10** until the `signed_circulation` origin/ordering contract (below) is
    nailed down and a throwaway spike confirms the correction actually fires and inverts the
    re-encoded marker on the real fixture.

---

## Section 3 (Scope) Analysis

* **Clarifying Questions:**
  - `OINVector` objects are only created inside the template branch (`oin_parser.py:495`).
    Ferrocene is `Fe_LIN` (2 slots) so vectors exist — but is there **any** in-scope eta
    geometry whose template lacks the slot, where `OINVector.winding` would be empty at the
    placement site even though winding is in-scope? If so the correction silently no-ops.
  - "Consume the single non-`None` `.winding` in the slot-group." What is the defined behavior
    when **zero** atoms in the group carry winding (a legacy V3.0–V3.5 OIN, or a hand-authored
    OIN with an eta ring but no marker)? The implied answer is "no correction," which means the
    generated face is arbitrary for those inputs — is that acceptable and documented?
  - What happens when **two** atoms in one slot-group carry conflicting non-`None` winding
    (malformed input)? "Single non-None value" assumes exactly one; the guard at
    `oin_parser.py:492` last-writer-wins on `winding_by_slot` but the per-vector `.winding`
    list retains both.

* **What-If Scenarios:**
  - **Malformed multi-marked ring.** Input `...c{0>}...c{0<}...` (two markers, same slot).
    `winding_by_slot` keeps the last; the vector list keeps both; the placement site picks
    whichever it iterates first. Encode and generate disagree on which marker "won" → silent
    face inversion with no error.
  - **Winding on a template-less in-scope path.** If any eta fixture routes through a
    geometry without a matching template slot, no `OINVector` is emitted, the correction is
    skipped, and US-001 passes for ferrocene while quietly failing for that shape.

* **Points for Improvement:**
  - Add a **negative/degenerate input contract** to scope: (a) zero-marker eta ring →
    documented no-op; (b) multi-marker same-slot → raise or canonicalize, do not silently pick.
  - Add an explicit in-scope assertion that **every in-scope eta geometry emits `OINVector`s
    carrying winding** (guard against the `oin_parser.py:495` template-gating hole), or move
    winding consumption to a channel that survives template-less paths.

---

## Section 4 (User Stories) Analysis

* **Clarifying Questions:**
  - **US-001** asserts a per-ring regex round-trip. Which atom does the re-encoder pick as the
    heading/star atom, and is that choice stable enough that the regex anchors to a *ring* not
    a *specific atom*? (R2 concedes this is the residual risk.) If the heading atom moves, the
    `{n>}`/`{n<}` may attach to a different ring atom and the regex must not care.
  - **US-002** parity test: "a placed geometry's helper-flip actually inverts the re-encoded
    marker." Does "helper-flip" mean applying the 180° rotation, or negating the helper's
    output? These test different things — the former validates the operator, the latter only
    validates arithmetic. Which is intended?
  - **US-005** requires "non-eta XYZ output byte-unchanged." `test_winding_inertness.py`
    already exists — does it already provide this harness, or is a new differ required? Is
    byte-identity robust to the numpy/scipy import reordering the refactor may introduce?

* **What-If Scenarios:**
  - **US-004 "conflict" is undefined.** "Both rings disagree in the same sense" vs "conflict"
    — what precisely is a conflict? Ring A wants a flip, Ring B does not. The PRD says leave
    as-is and log. But "no regression vs current output" is only true if *current* output is
    also unflipped; if today's code happens to place Ring A correctly by luck, the new
    "leave-as-is" path is a regression for the other ring. The baseline for "no regression"
    is unpinned.
  - **US-003 skip leaks to runtime.** The ferrocene *test* is skipped, but a user can still
    feed a symmetric-ring OIN with a winding marker at runtime. The correction will measure a
    geometrically-defined-but-physically-meaningless circulation on 5 identical carbons and may
    rotate pointlessly, then re-encode an arbitrary marker. US-003 covers the suite, not the
    runtime path.

* **Points for Improvement:**
  - Give **US-004 a concrete conflict definition and a pinned baseline artifact** (a committed
    "before" XYZ for a bridged ansa-metallocene) so "no regression" is falsifiable.
  - Add an acceptance criterion that the correction is a **no-op (identity, not a wasted
    rotation) on a detected-symmetric ring**, so runtime symmetric input is inert, matching the
    "not a geometric observable" resolution.
  - Split US-002 into two: (a) operator correctness (rotation inverts re-encoded marker on a
    real placed ring), (b) helper self-consistency (encode-side vs generate-side agree on a
    shared synthetic geometry).

---

## Section 5 (Technical Specifications) Analysis

* **Clarifying Questions:**
  - **`signed_circulation` origin contract (highest priority).** §5.4 signature is
    `signed_circulation(coords, star_local_idx, axis)` with `coords` = "placed binding-atom
    positions." The encoder centers first: `v_star = coord_star - centroid`
    (`oin_aligner.py:676`). Does the helper subtract the centroid of `coords` internally? If
    not, `_determine_winding` cannot be a faithful call site without pre-centering, and the two
    sides diverge — recreating the drift bug. **The centroid must be part of the helper's
    contract, or the API must take pre-centered vectors.** This is not documented.
  - **`v_next` ordering.** The encoder's "next" is the cyclic +1 in `constituent_indices`
    (sorted by `zone_a_info[k][3]`, `oin_aligner.py:661`). §5.4 says `coords` are in
    "fragment/SMILES order." Are these orders provably identical for eta rings? If the sort key
    reorders ring atoms relative to SMILES order, `v_next` picks a different neighbor and the
    sign flips. Prove or normalize.
  - **`star_local_idx` identity.** The encode side's star is the OIN heading atom; the generate
    side must locate the *same* atom's index within `coords`. What maps the parsed heading atom
    to the local index in the placed `binding_idxs`? If the placement reorders binding atoms
    (the analytic fallback reorders positions into `[C,H,C,H,…]` then re-stacks —
    `molassembler_adapter.py:723-732`), `star_local_idx` must track that reordering.
  - **Where does the `n < 3 → '>'` guard live** after the refactor? The helper returns only
    `'>'`/`'<'`. If the guard stays in `_determine_winding` but not in the generate caller, the
    two sides behave differently for degenerate rings.
  - **Axis definition parity.** Encoder dots with `slot_z` ("usually normalized or Pos vector,"
    `oin_aligner.py:691`); generate uses metal→centroid. The sign of `cross·axis` is invariant
    under proper rotation *and* positive scaling but **flips under axis negation**. Is
    metal→centroid guaranteed to point the same way as the encoder's `slot_z` (outward), for
    every slot, including slots where the template Z is defined inward?
  - **In-plane axis degeneracy.** The chosen axis is "centroid→binding-atom[0] projected into
    the ring plane." If binding-atom[0] sits near the metal→centroid axis (its in-plane
    projection ≈ 0), the rotation axis is ill-conditioned. What is the numerical guard and the
    fallback axis?

* **What-If Scenarios:**
  - **Silent sign drift via un-centered helper.** Implementer writes `signed_circulation` using
    raw `coords` (no centering) because the API doesn't mention a centroid. Encode still
    centers (existing code), generate uses raw. Both call "the one helper," the grep in §8
    passes, and the correction inverts the wrong inputs — the exact failure the shared helper
    was created to prevent, now *harder* to spot because it looks unified.
  - **Reorder-induced star mismatch.** Analytic fallback path fires (ETKDG fails to embed the
    substituted Cp anion). Positions get re-stacked to `[heavy…, H…]`
    (`:723-732`). `star_local_idx` computed against the original `binding_idxs` now points at
    the wrong row of `coords`. Correction rotates about a wrong axis / measures wrong sign.
  - **Axis-negation flip on a lower-hemisphere slot.** For an octahedral/`LIN` slot whose
    template Z points toward the metal, `slot_z` and metal→centroid are antiparallel; the
    helper returns opposite chars on the two sides for identical geometry.
  - **Non-planar substituted ring.** Heavy halides (I, Br) pucker the ring; the SVD plane
    normal and the true circulation axis diverge by several degrees. `cross(v_star,v_next)` is
    no longer parallel to `axis`; sign is still robust *unless* the ring is near-perpendicular
    to a bad axis choice — combine with the degeneracy case above and the sign becomes unstable.

* **Points for Improvement:**
  - **Rewrite the §5.4 contract** to pin: (a) `coords` are centered on their own centroid
    *inside* the helper (or the docstring states they must be pre-centered and both call sites
    do so identically); (b) the ordering of `coords` is exactly the SMILES/fragment order and
    `_determine_winding` is refactored to feed that same order; (c) the `n<3` default lives
    *inside* the helper so both sides inherit it; (d) `axis` orientation convention is fixed to
    metal→centroid (outward) and the encoder's `slot_z` feed is asserted to match, with a unit
    test on an antiparallel slot.
  - Add a **helper property test independent of the fixture**: synthesize a known ring
    (regular pentagon with one labeled substituent), feed identical coords to both call paths,
    assert identical char; then apply the 180° rotation and assert the char flips. This closes
    the drift class by construction rather than by grep.
  - Add a **guard + deterministic fallback** for the in-plane axis (e.g. if
    `‖projection‖ < ε`, use binding-atom[1], then Gram–Schmidt against `axis`).
  - Specify **idempotency**: applying the correction twice (or round-tripping
    generate→encode→generate) must converge, not oscillate. Add an assertion.

---

## Section 6 (Negative Constraints) Analysis

* **Clarifying Questions:**
  - "DO NOT duplicate the winding-sign math — call `signed_circulation`." The §8 success metric
    enforces this by **grep for cross/dot in `molassembler_adapter.py`**. But an un-centered or
    mis-ordered call to the shared helper passes grep while still being wrong. Is there a
    *behavioral* check, not just a textual one?
  - "DO NOT use a reflection operator — it inverts pendant chirality." The chosen fixture
    (H/OH/Cl/Br/I substituents) contains **no stereocenter**, so a mistakenly-reflecting
    implementation would pass every acceptance test. How is the `det +1` property actually
    verified?

* **What-If Scenarios:**
  - **Reflection ships undetected.** An implementer uses `Rotation.align_vectors` in a way that
    produces an improper transform, or negates a coordinate to "flip the face." All acceptance
    tests are green because the fixture has no chirality to invert. The bug detonates in Phase 4
    when a real ligand stereocenter is present.

* **Points for Improvement:**
  - Add a **`det(R) ≈ +1` runtime assertion** on the correction transform (cheap, exact,
    catches every improper-operator regression at the source).
  - Add a minimal **chirality-witness fixture** (one eta ring bearing a single pendant
    stereocenter) whose CIP code must be invariant across the correction — even if full Phase 4
    encoding is out of scope, a "chirality preserved" assertion is in-scope for *this*
    operator's correctness.

---

## Section 7 (Risks & Mitigation) Analysis

* **Clarifying Questions:**
  - R1 (sign inversion) is mitigated by "shared helper + parity test." Given the un-centered /
    ordering ambiguities in §5, does the parity test feed *raw placed coordinates from a real
    generation run* through both call sites, or a synthetic ideal? Only the former catches the
    reorder/centroid bugs.
  - R2 lists a "documented weaker fallback (differ + geometric halide-sequence-reversal
    check)." Who decides when the primary exact assertion is abandoned for the fallback — is the
    fallback a *test-time* human decision or an automated `assert-or-skip`? An auto-skip that
    silently downgrades hides real regressions.
  - R4 flags a 0.95 Å O–H and a "poorly-defined Kabsch warning" in the fixture. Has the fixture
    been re-verified to encode **deterministically** across runs, given the golden does not yet
    exist on disk?

* **What-If Scenarios:**
  - **Golden blessed from a wrong artifact.** `Ferrocene-halide-face_oin.txt` is created by
    running the not-yet-correct pipeline once and committing whatever it emits. If the
    correction has a sign bug, the golden enshrines the bug and the round-trip is "faithful to
    the wrong answer."
  - **Fallback masks R1.** The R2 fallback auto-triggers because the exact heading-atom
    assertion is flaky; it checks only the halide-sequence order, which a *reflection* would
    ALSO satisfy (sequence reverses under both proper-rotation-flip and reflection). R1 and the
    reflection ban both slip through the weakened check.

* **Points for Improvement:**
  - The golden must be **derived from an independently-reasoned expected string**, not
    machine-emitted-then-blessed. The PRD already asserts the expected encoding in §5.1 — pin
    the golden to *that* hand-verified string and fail if generation disagrees.
  - Make R2's fallback an **explicit `xfail`/documented-skip with a logged reason**, never a
    silent auto-downgrade, so a flaky primary assertion is visible in CI, not swallowed.
  - Add a risk **R6 — observability**: correction decisions currently surface only via
    `print(..., file=sys.stderr)` (the module's existing pattern). A caller cannot tell whether
    a face was corrected, skipped, or hit the multi-eta conflict path. Add a structured,
    inspectable signal (return metadata or a proper logger) so silent mis-placement (US-004
    conflict, zero-marker no-op) is detectable in production.

---

## Section 8 (Success Metrics) Analysis

* **Clarifying Questions:**
  - The grep metric ("no duplicated cross/dot sign logic") is a **textual** proxy for a
    **behavioral** property. What behavioral metric proves encode and generate actually agree
    on a shared placed geometry?
  - "Full `unittest discover tests/unit` green" — does the metric include the pre-existing
    `test_winding_inertness.py`, and is that suite currently green on `main` so a regression is
    attributable to this change?
  - "Non-eta XYZ output byte-unchanged" — byte-identity or geometric identity? The refactor
    adds imports and a parameter; if any non-eta path formatting shifts, byte-diff fails
    spuriously. Is there tolerance, or is strict byte-identity the contract?

* **What-If Scenarios:**
  - **Green suite, wrong physics.** Every listed metric passes with an un-centered helper and a
    stereocenter-free fixture (per §5/§6 findings). The metrics as written do not distinguish
    "correct" from "self-consistently wrong."

* **Points for Improvement:**
  - Add these metrics: (1) **behavioral parity** — one placed ring's helper char equals the
    encoder's char on the identical coordinate array; (2) **operator invert** — applying the
    correction flips the re-encoded char (instrumented, both branches); (3) **`det +1`** on the
    correction transform; (4) **CIP invariance** on a chirality-witness ring; (5) **idempotency**
    on a double-applied correction; (6) **golden equals the §5.1 hand-verified string**, not
    merely "round-trips."

---

## Consolidated Top Findings (ranked for `/hyper-resolve`)

1. **`signed_circulation` centroid/ordering/axis contract is underspecified** (§5) — the single
   highest risk; an un-centered or mis-ordered "shared" helper reintroduces the sign-drift class
   while passing the §8 grep. *Fix: pin origin, order, `n<3` default, and axis orientation in the
   API contract; add a fixture-independent behavioral parity test.*
2. **Chosen fixture cannot test the `det +1` justification** (§6) — no pendant stereocenter, so a
   reflection bug ships green and detonates in Phase 4. *Fix: `det(R)≈+1` assertion + a
   chirality-witness ring with CIP invariance.*
3. **Fixed ETKDG seed may leave one correction branch uncovered** (§1) — flip test may only
   exercise the rotation in one direction. *Fix: instrumented two-branch coverage.*
4. **Golden does not yet exist and risks self-blessing a wrong artifact** (§7) — *Fix: pin golden
   to the §5.1 hand-verified string.*
5. **US-004 "conflict" undefined and "no regression" baseline unpinned** (§4) — *Fix: concrete
   conflict definition + committed baseline XYZ.*
6. **No observability of correction decisions** (§7, R6) — silent no-op / conflict paths are
   invisible. *Fix: structured signal instead of stderr prints.*
7. **Degenerate in-plane axis has no specified guard** (§5) — *Fix: ε-guard + deterministic
   fallback axis.*
8. **Winding channel gated on template branch** (§3, `oin_parser.py:495`) — confirm no in-scope
   eta geometry loses `OINVector.winding`. *Fix: assert emission or move channel.*

**Blast radius (confirmed contained):** `_stitch_eta_fragment` and `_determine_winding` each have
a single production caller; the new parameter defaults to `None` (back-compat); the `:1170`
optimiser is rotation-only and preserves circulation. The refactor's principal external hazard is
the pre-existing `test_winding_inertness.py` inertness contract (US-005), which must stay green.

**Recommended confidence after triage: 7/10** (from the draft's 9/10) until findings 1, 2, and 4
are resolved.

---

**Next step:** Run `/hyper-resolve` to triage these findings and compile the final SuperPRD +
MiniPRDs (`MiniPRD_SignedCirculationHelper_Phase3.md`, `MiniPRD_HapticFaceCorrection_Phase3.md`).
