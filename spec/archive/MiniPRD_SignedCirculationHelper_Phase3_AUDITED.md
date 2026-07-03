# MiniPRD: SignedCirculation Helper (Stereo Phase 3)
**Hypergraph Node ID:** oin.winding.signed_circulation
**Parent Node:** SuperPRD_StereoPhase3_HapticFace (Stereo Phase 3 — Haptic Face Correction)

## 1. The Confidence Mandate
**Confidence Score**: 10/10.

Mechanical extraction of an existing, working sign computation into a shared, fully-contracted
helper. The one prior ambiguity (centroid/ordering/axis) is now pinned by the Red Team resolution;
no residual clarifying questions.

## 2. Atomic User Stories
* **US-002 (from SuperPRD):** As a maintainer, I want encode and generate to share one
  winding-sign definition, so the correction cannot silently invert due to convention drift.
  * The helper is the single source of truth; both `_determine_winding` (encode) and the
    generation-side correction call it.
  * Behavioral parity: a synthetic placed ring's helper char equals the encoder's char on the
    identical coordinate array (fixture-independent).

## 3. Implementation Plan (Task List)
- [ ] Create `src/oinsmiles/oin/winding.py` with `signed_circulation(coords, star_local_idx, axis) -> str`.
- [ ] Inside the helper: center `coords` on their own centroid (`coords - coords.mean(0)`) — **do not** assume pre-centered input.
- [ ] Inside the helper: host the `n < 3 → '>'` default so both call sites inherit it.
- [ ] Inside the helper: `v_star = centered[star_local_idx]`; `v_next` = the cyclic +1 neighbour in the given (SMILES/fragment) order; return `'>'` if `cross(v_star, v_next)·axis >= 0` else `'<'`.
- [ ] Document in the docstring: coords MUST be SMILES/fragment order; axis convention is metal→centroid **outward**.
- [ ] Refactor `_determine_winding` (`oin_aligner.py:619`, sole caller) to build coords in the **same SMILES/fragment order** and delegate the sign to `signed_circulation` — remove the inline `cross`/`dot` sign math (`oin_aligner.py:676-691`).
- [ ] Assert/normalize that the encoder's `slot_z` feed points outward (metal→centroid), matching the helper's axis convention.
- [ ] Add an **antiparallel-slot unit test**: an octahedral/`LIN` slot whose template Z points inward must still yield the correct char (guards axis-negation flip).
- [ ] Add a **fixture-independent behavioral parity test**: synthesize a regular pentagon with one labeled substituent; feed identical coords to both the encode path and a direct helper call → assert identical char; then apply the 180° in-plane rotation → assert the char flips.

## 4. The Negative Space (Constraints)
* **DO NOT** accept raw un-centered coords or a non-SMILES order — centering and ordering are part of the helper's contract.
* **DO NOT** leave the `n < 3` default in `_determine_winding` only — it must live inside the helper.
* **DO NOT** leave duplicated `cross`/`dot` sign logic in `oin_aligner.py` or `molassembler_adapter.py` after the refactor (grep must be clean).
* **DO NOT** silently rely on `slot_z` being outward — assert the convention.

## 5. Integration Tests & Verification
* **Test 1 (Deterministic — parity):** identical synthetic pentagon coords → encode-path char == direct `signed_circulation` char. Expected: equal.
* **Test 2 (Deterministic — operator invert):** apply the 180° in-plane rotation to the synthetic ring → `signed_circulation` returns the opposite char. Expected: flipped.
* **Test 3 (Deterministic — antiparallel slot):** inward-Z slot geometry → char matches the outward-convention expectation. Expected: convention-correct.
* **Test 4 (Deterministic — degenerate):** `n < 3` ring → `'>'`. Expected: `'>'`.
