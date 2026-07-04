# SuperPRD — Eta-Ring Canonicalization (Fragment Order + Heading Atom)

## Metadata
- **Project Name**: OIN-SMILES — Substituted Eta-Ring Canonical-Form Stability (Stereo Phase 3 residual)
- **Version**: 1.0.0
- **Status**: Approved (Red-Team resolved)
- **Owner**: Resolution Agent (HACF chain: architect → redteam → resolve)
- **Target Node**: `atom_oin_aligner` (`src/oinsmiles/utils/oin_aligner.py`)
- **Confidence Score**: 10/10 (both blocking findings RT-1, RT-4 resolved with explicit decisions)

## 1. Introduction & Goals

### 1.1 Problem Statement
The XYZ→OIN encoder's canonical form is **not invariant** across two physically-equivalent
inputs of a substituted-eta-ring complex (a hand-built fixture vs. an ETKDG-generated
structure of the same molecule). `Ferrocene-halide-face.xyz` (ferrocene whose two Cp rings
carry *different* pentahalo substitutions, `O,H,Cl,Br,I` and `O,H,I,Br,Cl`) round-trips
XYZ→OIN→generate-3D→re-encode to a **different byte string** than its pinned golden, though
the physical winding is preserved. This is the sole remaining `@unittest.expectedFailure`:
`tests/unit/test_stereo_roundtrip_diagnostics.py::test_haptic_face_golden_match`.

Measured diff (real data — do not re-derive):
```
GOLDEN  : [Fe_LIN].Oc{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I . Oc{1}1[cH]{1}c{1}(I)c{1<}(Br)c{1}1Cl
REENCODE: [Fe_LIN].Oc{0}1[cH]{0<}c{0}(I)c{0}(Br)c{0}1Cl . Oc{1}1[cH]{1}c{1}(Cl)c{1}(Br)c{1<}1I
```
Matched by CONTENT, this is two winding-preserving relabelings: **RC1** the two distinct
rings swap fragment order; **RC2** the winding marker drifts to a different ring atom within a
ring of fixed content (same traversal direction, same marker CHARACTER).

### 1.2 Solution Overview
Make the canonical form invariant to (RC1) eta-ring fragment arrival order and (RC2) geometric
ring orientation for substituted eta rings, as a **pure winding-preserving relabeling** — only
WHICH ring is first and WHICH atom is the visible heading change; the winding CHARACTER
(`>`/`<`) is still computed from geometry by `_determine_winding`/`signed_circulation`
(FROZEN). Both fixes are internal to `OINDiscreteAligner`.

- **RC1 (fragment order):** a **scoped eta-only rank swap** — permute only same-mass eta
  fragments among the rank slots they already occupy, keyed on a **heading-independent RDKit
  canonical ring SMILES** with **winding sense as a content-identical tiebreak only**. Every
  non-eta fragment keeps its exact rank; the metal stays rank 0. (D-RC1; RT-1 disposition.)
- **RC2 (heading atom):** for substituted (asymmetric) eta rings, heading = the ring atom with
  the **lowest `Chem.CanonicalRankAtoms(..., breakTies=True)`** — topological, order-invariant,
  the content-derived analog of the existing `SYMMETRIC_LIGANDS` "lowest local index" override
  (which stays first-wins for symmetric rings). Guarded by a **proven start-invariance
  property** (RT-4 disposition). (D-RC2.)

### 1.3 Target Audience
Internal — the XYZ→OIN encoder and its golden-string contract. No public API change.

## 2. Confidence Mandate
**Confidence Score**: 10/10. Root causes pinned to file:lines; both design decisions resolved;
both blocking Red-Team findings (RT-1, RT-4) dispositioned with concrete engineering
directives; the byte-identical golden set is the objective acceptance gate.

**Clarifying Questions**: all resolved (see §5.1 Resolved Trade-offs Log).

## 3. Scope

### 3.1 In-Scope
- RC1 scoped eta-only rank swap for same-mass eta fragments.
- RC2 canonical-rank heading for substituted/asymmetric eta rings.
- Un-`expectedFailure` and pass `test_haptic_face_golden_match`.
- New hard tests per RT-5 (inertness, scoped-swap, start-invariance, live reflection).

### 3.2 Out-of-Scope
- Any change to the winding-sign computation (`_determine_winding` / `signed_circulation`) — FROZEN.
- The `SYMMETRIC_LIGANDS` forced-heading path; all non-eta / heterogeneous-eta fragment order.
- Generation (OIN→XYZ) side, molassembler adapter, Zone-A P/N stereo.
- Deleting/weakening/auto-wiring the R2 geometric-fallback skip.

## 4. User Stories (Atomic)
| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | Substituted-eta-ring FRAGMENTS ordered by content-canonical key, independent of xyz2mol arrival order. | 1. Two content-distinct eta rings order identically for fixture and generated inputs.<br>2. Key = heading-independent canonical ring SMILES; winding tiebreak only for content-identical rings.<br>3. ONLY same-mass eta fragments reorder; only within the rank slots they already occupy; every non-eta rank byte-identical; metal stays rank 0. | High |
| US-002 | Heading atom of a substituted eta ring chosen from ring TOPOLOGY (lowest canonical rank), not 3D orientation. | 1. Heading identical for fixture and generated structures of the same ring.<br>2. `SYMMETRIC_LIGANDS` path untouched (first-wins).<br>3. Winding character still from `_determine_winding`/`signed_circulation`. | High |
| US-003 | The golden round-trip closes. | 1. `test_haptic_face_golden_match` un-`expectedFailure`d and PASSES byte-for-byte.<br>2. `discover tests/unit` → skipped=3, **expected failures=0**.<br>3. Every pre-existing golden byte-identical. | High |
| US-004 | Canonicalization never masks a real reflection. | 1. Start-invariance proven: character identical for every star choice on the fixture ring.<br>2. A reflected ring still yields a winding-flipped marker (live assertion).<br>3. R2 skip stays present, skipped, and meaningful. | Critical |

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Resolved Trade-offs

**Data flow (unchanged):** `xyz2mol.get_tmc_mol()` → `CIPAssigner.assign_all()` →
`OINDiscreteAligner._reduce_hapticity()` (builds virtual-atom dicts, sets `"rank": i`) →
`_permute_and_serialize()` (homogeneous sort + heading selection + winding) → V3.6 string.

**Root-cause mechanism (corrected by RT-1):**
- The serialized fragment order comes from `final_sorted_view.sort((rank, local_idx))`
  (`:497`), where `rank` originates as `"rank": i` in `_reduce_hapticity` (`:306`). The
  homogeneous sort (`:456-494`) only permutes items *within an equal `chem_id` bucket*, so two
  **content-distinct** rings (different `chem_id`) are never compared — their order is the
  arrival index `i`. `base_sort_key` (`:250`, stored as `"key"` at `:308`) is **dead code**
  (never read) and is NOT the RC1 lever.
- RC2 heading (`:508-563`) is a geometric `best_idx` (max alignment of centroid→atom to the
  slot `ref_vec`); the `SYMMETRIC_LIGANDS` override (`:565-597`) forces lowest local index but
  does not fire for substituted rings.

**Resolved Trade-offs Log:**

- **RT-1 (BLOCKING) — RC1 lever & golden-shift risk.** *Options:* (A) scoped eta-only rank
  swap; (B) global re-rank with golden-preserving tiebreak. *Resolution:* **(A).** Permute only
  same-mass eta fragments among the rank slots they already occupy; all non-eta fragments and
  the metal (rank 0) keep their exact rank. A global re-rank was rejected: canonical `N` sorts
  before `[Cl]`, so it inverts cisplatin/transplatin — unacceptable golden shift. The dead
  `base_sort_key` is retired; RC1 operates on the rank-assignment step, not `chem_id`.
- **RT-4 (BLOCKING) — winding start-invariance / reflection masking.** *Options:* (A)
  prove-and-assert, keep winding math frozen; (B) rewrite the character to be explicitly
  heading-independent. *Resolution:* **(A).** `signed_circulation`/`_determine_winding` stay
  FROZEN. For the substituted-eta construction, `constituent_indices = sorted(local_idx)` of
  the ring atoms tracks ring-cyclic order and the ring is planar-convex, so the single-edge
  sign is start-invariant. This is made **load-bearing via tests**: (i) a hard assertion that
  the character is identical for every choice of star on the fixture ring, and (ii) a live
  reflection test that a mirrored ring still flips the marker. Option B was rejected: it
  unfreezes the winding math and touches `signed_circulation`'s shared contract with the
  generation-side correction. If start-invariance ever fails to hold, the fix does not land.
- **D-RC2 (heading rule).** Lowest `Chem.CanonicalRankAtoms(breakTies=True)`; topological,
  order-invariant, smallest delta from the existing symmetric override.
- **D-RC1 (order identity).** Heading-independent canonical ring SMILES primary; winding sense
  tiebreak confined to content-identical rings (for the halide fixture the rings differ in
  content, so winding never enters the key — safety preserved).
- **RT-2 (HAZARD) — no mol at the site.** *Resolution:* build/thread a fragment mol whose atom
  indices map back to `constituent_indices`/`local_idx` without a round-trip; prefer an
  atom-map-preserving construction over bare `MolFromSmiles(lig["smiles"])`. Compute canonical
  rank/SMILES from that mol.
- **RT-3 (HAZARD) — unstable `lig["smiles"]` & bond-perception divergence.** *Resolution:* the
  RC1 signature must be a fresh `CanonicalRankAtoms`/`MolToSmiles` (order-invariant), computed
  after the same perception path that produces `lig["smiles"]` so fixture and generated
  fragments perceive bonds identically; on any failure, fall back to today's behavior for that
  fragment (fail-safe, never a silent partial reorder).
- **RT-5 (HARDEN) — missing tests + incomplete key.** *Resolution:* add the four hard tests
  (§8); for content-identical + same-winding rings (key still ties), add a final deterministic
  tiebreak (lowest constituent global index) so no arrival-order leak remains.

### 5.2 System Graph Blast Radius
- **Modified:** `atom_oin_aligner` (`src/oinsmiles/utils/oin_aligner.py`) — ONLY node.
- **Must stay byte-inert:** `atom_oin_writer`, `atom_oin_sanitizer`, `atom_cip_assigner`,
  `atom_xyz2mol`, `oin.winding.signed_circulation` (FROZEN), entire generation side.

### 5.3 Execution Checklist (MiniPRDs)
- [ ] `spec/compiled/MiniPRD_EtaRingCanonicalization.md`

### 5.4 API Contracts / Schema
No public API change. Internal to `OINDiscreteAligner`:
- RC1: a scoped rank-permutation over same-mass eta fragments, keyed
  `(canonical_ring_smiles, winding_sense_tiebreak, lowest_constituent_global_idx)`; permutes
  only within the eta fragments' existing rank-slot set.
- RC2: heading = `min` by `Chem.CanonicalRankAtoms(mol, breakTies=True)` over the ring's
  constituent atoms, gated behind the `SYMMETRIC_LIGANDS` first-wins check; falls back to the
  current geometric `best_idx` if the canonical rank cannot be computed/mapped.

### 5.5 Dependencies
- RDKit (`Chem.CanonicalRankAtoms`, `Chem.MolToSmiles`) — existing hard dependency. No new libs.

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** modify `_determine_winding` or `oin.winding.signed_circulation` (FROZEN).
- **DO NOT** touch the `SYMMETRIC_LIGANDS` forced-heading path or any non-eta / heterogeneous
  fragment ordering — gate all new logic on "substituted/asymmetric eta group."
- **DO NOT** re-rank any non-eta fragment or move the metal off rank 0.
- **DO NOT** let winding sense enter the RC1 key except as the content-identical tiebreak.
- **DO NOT** delete, weaken, un-skip, or auto-substitute
  `test_haptic_face_r2_geometric_fallback_never_auto_substituted`.
- **DO NOT** re-pin any existing golden; the fix must reproduce them byte-for-byte.

## 7. Risks & Mitigation
- **R1 — golden shift from shared path** → scoped eta-only rank swap (RT-1/A); full golden set
  is a hard gate.
- **R2 — reflection masking (CRITICAL)** → winding math frozen; start-invariance proven by
  test; winding in RC1 key only as content-identical tiebreak; R2 skip stays; live reflection
  test added (RT-4/A).
- **R3 — mol/index provenance & bond-perception divergence** → correct mol source + round-trip-
  free index map; same-perception signature; fail-safe fallback (RT-2/RT-3).

## 8. Success Metrics
- `test_haptic_face_golden_match` un-`expectedFailure`d and PASSES byte-for-byte against
  `tests/candidate_outputs/Ferrocene-halide-face_oin.txt`.
- `uv run python -m unittest discover tests/unit` → OK, **skipped=3, expected failures=0**.
- `uv run python -m unittest discover tests` → OK.
- New hard tests pass: (1) non-eta order inertness (cisplatin/transplatin/BDPP/BDNN/BINAP);
  (2) symmetric-eta inertness (plain ferrocene / ansa-metallocene); (3) RC1 scoped-swap
  (eta ring at a non-1 rank keeps non-eta ranks fixed); (4) RC2 start-invariance + a reflected
  ring still flips the marker.
- Per-ring content-anchored tests still pass; R2 skip present and meaningful.

## 9. Red-Team Disposition Log
| ID | Severity | Disposition |
|----|----------|-------------|
| RT-1 | BLOCKING | Scoped eta-only rank swap (Option A). base_sort_key retired; fix at rank-assignment; non-eta ranks + metal untouched. |
| RT-4 | BLOCKING | Prove-and-assert; winding math FROZEN (Option A). Start-invariance asserted by test; live reflection test added. |
| RT-2 | HAZARD | Correct mol source + round-trip-free index map to `constituent_indices`. |
| RT-3 | HAZARD | Order-invariant same-perception signature; fail-safe fallback on any compute failure. |
| RT-5 | HARDEN | Four hard tests added; final deterministic tiebreak (lowest constituent global idx) for identical-content/same-winding rings. |
