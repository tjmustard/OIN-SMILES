# Draft PRD — Eta-Ring Canonicalization (Fragment Order + Heading Atom)

## Metadata
- **Project Name**: OIN-SMILES — Substituted Eta-Ring Canonical-Form Stability (Stereo Phase 3 residual)
- **Version**: 0.1.0 (Draft)
- **Status**: Draft
- **Owner**: Architect Agent (seeded from `spec/worklog/haptic-canonicalization-design-brief.md`)
- **Target Node**: `atom_oin_aligner` (`src/oinsmiles/utils/oin_aligner.py`)

## 1. Introduction & Goals

### 1.1 Problem Statement
The XYZ→OIN encoder produces a canonical string that is **not invariant** across two
physically-equivalent inputs of a substituted-eta-ring complex (a hand-built fixture vs. an
ETKDG-generated structure of the same molecule). Concretely, `Ferrocene-halide-face.xyz`
(a ferrocene whose two Cp rings carry *different* pentahalo substitutions:
`O,H,Cl,Br,I` and `O,H,I,Br,Cl`) round-trips XYZ→OIN→generate-3D→re-encode to a
**different byte string** than its pinned golden, even though the physical winding is
preserved. This is the sole remaining `@unittest.expectedFailure` in the suite:
`tests/unit/test_stereo_roundtrip_diagnostics.py::test_haptic_face_golden_match`.

Measured diff (2026-07-03, real data — do not re-derive):
```
GOLDEN  : [Fe_LIN].Oc{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I . Oc{1}1[cH]{1}c{1}(I)c{1<}(Br)c{1}1Cl
REENCODE: [Fe_LIN].Oc{0}1[cH]{0<}c{0}(I)c{0}(Br)c{0}1Cl . Oc{1}1[cH]{1}c{1}(Cl)c{1}(Br)c{1<}1I
```
Matching rings by CONTENT (not position) shows exactly two winding-preserving relabelings,
neither reversing winding:
1. **RC1 — fragment-order swap.** The two distinct rings appear in opposite fragment order.
2. **RC2 — heading-atom drift.** Within a ring of fixed content, the winding marker `<` sits
   on a different ring atom between the two encodings; traversal direction and marker
   CHARACTER are identical → winding sense unchanged, only the visible heading atom moved.

### 1.2 Solution Overview
Make the encoder's canonical form invariant to (RC1) eta-ring fragment arrival order and
(RC2) geometric ring orientation for substituted eta rings, as a **pure
winding-preserving relabeling** — the fix changes only WHICH ring is listed first and WHICH
atom is the visible heading; it NEVER changes how the winding character (`>`/`<`) is computed.
Both fixes live in `OINDiscreteAligner` (`src/oinsmiles/utils/oin_aligner.py`):

- **RC1 (fragment order):** replace the arrival-order index `i` in the fragment sort identity
  (`base_sort_key`/`chem_id`, `:250`, `:313`) with a **heading-independent RDKit canonical
  ring SMILES**; use **winding sense only as a final tiebreak** between content-identical
  rings. (Decision D-RC1, §5.1.)
- **RC2 (heading atom):** for substituted (asymmetric) eta rings, replace the
  geometry-dependent `best_idx` (`:540-560`) with a content-canonical heading = the ring atom
  with the **lowest RDKit canonical atom rank** (`Chem.CanonicalRankAtoms`), analogous to the
  existing `SYMMETRIC_LIGANDS` forced-heading override (`:565-597`, "lowest local index") but
  derived from topology instead of arrival index. (Decision D-RC2, §5.1.)

### 1.3 Target Audience
Internal — the XYZ→OIN encoder pipeline and its golden-string contract. No public API change.

## 2. Confidence Mandate
**Confidence Score**: 8/10
The root causes are pinned to exact file:lines and both design decisions are resolved
(§5.1). Residual uncertainty is implementation-local: (a) whether the RDKit canonical ring
SMILES/rank is computable from the exact mol object available at the RC1/RC2 sites without a
re-sanitize round-trip, and (b) confirming no golden other than the target shifts (the shared
`base_sort_key` is the blast-radius concern — see §7 R1). Both are settled at execute time by
the byte-identical golden gate, not by more human input.

**Clarifying Questions** (both RESOLVED in this session — see §5.1):
- [x] RC2 heading rule → RDKit canonical atom rank (lowest = heading).
- [x] RC1 fragment-order identity → canonical ring SMILES, winding sense as tiebreak only.

## 3. Scope

### 3.1 In-Scope
- RC1: content-canonical fragment-order identity for **same-mass eta groups only**
  (gate any new logic on "the fragment is an eta/haptic group"; heterogeneous or non-eta
  fragments unchanged).
- RC2: content-canonical heading-atom selection for **substituted/asymmetric eta rings only**
  (the `SYMMETRIC_LIGANDS` forced-heading path stays exactly as-is).
- Flipping `test_haptic_face_golden_match` from `@unittest.expectedFailure` to a real pass.

### 3.2 Out-of-Scope
- Any change to the winding-sign computation (`_determine_winding` / `signed_circulation`).
  The winding character MUST continue to be derived from geometry.
- The `SYMMETRIC_LIGANDS` forced-heading path and all non-eta / heterogeneous-eta fragments.
- The generation (OIN→XYZ) side, molassembler adapter, and Zone-A P/N stereo.
- Deleting, weakening, or auto-wiring the R2 geometric-fallback skip
  (`test_haptic_face_r2_geometric_fallback_never_auto_substituted`).

## 4. User Stories (Atomic)
| ID | User Story | Acceptance Criteria | Priority |
| :--- | :--- | :--- | :--- |
| US-001 | As the encoder, I want substituted-eta-ring FRAGMENTS ordered by a content-canonical key so fragment order is independent of xyz2mol arrival order. | 1. Two content-distinct eta rings order identically for fixture and generated inputs.<br>2. Order key = heading-independent canonical ring SMILES; winding sense only tiebreaks content-identical rings.<br>3. Only same-mass eta fragments are reordered; every other fragment order is byte-identical. | High |
| US-002 | As the encoder, I want the heading atom of a substituted eta ring chosen from ring TOPOLOGY (lowest RDKit canonical rank), not 3D orientation. | 1. Heading atom is identical for fixture and generated structures of the same ring.<br>2. `SYMMETRIC_LIGANDS` path untouched.<br>3. Winding character still computed by `_determine_winding`/`signed_circulation`. | High |
| US-003 | As a maintainer, I want the golden round-trip to close. | 1. `test_haptic_face_golden_match` un-`expectedFailure`d and PASSES byte-for-byte.<br>2. `discover tests/unit` → skipped=3, **expected failures=0**.<br>3. Every pre-existing golden byte-identical. | High |
| US-004 | As a safety reviewer, I want the canonicalization to never mask a real reflection. | 1. A genuinely reflected ring face still yields a winding-flipped (`>`↔`<`) marker.<br>2. The R2 skip stays present, skipped, and meaningful (not auto-substituted for the exact-match assertion). | Critical |

## 5. Technical Specifications (The Blueprint)

### 5.1 Architecture & Resolved Trade-offs

**Data flow (unchanged):** `xyz2mol.get_tmc_mol()` → `CIPAssigner.assign_all()` →
`OINDiscreteAligner._permute_and_serialize()` (canonical string) → V3.6 inline string. The
fix touches only the sort-identity construction and the heading-atom selection inside
`OINDiscreteAligner`.

**Current mechanism (root causes, both `oin_aligner.py`):**
- RC1: `_reduce_hapticity` builds `base_sort_key = (i, first_binding_atom_mass, lig["smiles"])`
  (`:250`) and `chem_id = (first_binding_atom_mass, lig["smiles"])` (`:313`). Homogeneous
  sorting (`_permute_and_serialize`, `:456-494`) groups by `chem_id` and reorders only WITHIN
  a `chem_id` group. Two content-distinct rings land in different `chem_id` groups, so they
  are never compared — their final `rank` is the arrival index `i`, which differs between
  fixture and generated. `lig["smiles"]` is itself heading/traversal-dependent, so it cannot
  serve as a stable identity as-is.
- RC2: heading selection (`:508-563`) picks `best_idx` = the ring atom whose
  centroid→atom vector, rotated into the template frame, maximally aligns with the slot
  `ref_vec` — a GEOMETRIC choice keyed on absolute 3D orientation. The `SYMMETRIC_LIGANDS`
  override (`:565-597`) forces heading = lowest constituent index but does NOT fire for
  substituted rings, which fall through to the unstable geometric pick.

**Resolved Trade-offs Log:**

- **D-RC2 (heading rule).** *Issue:* substituted eta rings need a heading atom that is
  identical from fixture and generated topology. *Options:* (A) RDKit canonical atom rank;
  (B) highest-CIP-priority substituent carbon; (C) canonically-first substituent SMILES.
  *Resolution:* **(A) lowest RDKit canonical atom rank** (`Chem.CanonicalRankAtoms`). Purely
  topological, already available in the pipeline, and is the content-derived analog of the
  existing `SYMMETRIC_LIGANDS` "lowest local index" rule — smallest conceptual delta. CIP
  coupling (B) was rejected as unnecessary coupling to the CIP pipeline; (C) needs a secondary
  key for identical substituents that (A) gets for free.

- **D-RC1 (fragment-order identity).** *Issue:* eta-ring fragments must order by content, not
  arrival. *Options:* (A) canonical ring SMILES with winding sense as tiebreak; (B) canonical
  ring SMILES only; (C) sorted substituent multiset + winding. *Resolution:* **(A)
  heading-independent RDKit canonical ring SMILES as primary key, winding sense used ONLY as a
  final tiebreak between content-identical rings.** Canonical SMILES is a stronger identity
  than a substituent multiset (C). Winding-in-key is confined to the content-identical case:
  for the halide fixture the two rings differ in content so winding never enters the key,
  preserving the safety property (§7 R2); (B) would tie two content-identical enantiomeric
  faces and need an extra fallback anyway.

- **Winding computation is FROZEN.** Both fixes only choose WHICH ring is first and WHICH atom
  is heading. The `>`/`<` character is still produced by `_determine_winding` →
  `signed_circulation` from geometry. This is the load-bearing safety invariant (US-004).

### 5.2 System Graph Blast Radius
The following nodes in `spec/compiled/architecture.yml` are affected:
- `atom_oin_aligner` (`src/oinsmiles/utils/oin_aligner.py`) — the ONLY node modified.

Untouched by construction (gate every new branch on "substituted/asymmetric eta group"):
`atom_oin_sanitizer`, `atom_cip_assigner`, `atom_xyz2mol`, `atom_oin_writer`,
`oin.winding.signed_circulation`, and the entire generation side.

### 5.3 Execution Checklist (MiniPRDs)
- [ ] `spec/compiled/MiniPRD_EtaRingCanonicalization.md`

### 5.4 API Contracts / Schema
No public API change. Internal, within `OINDiscreteAligner`:
- The eta-fragment sort identity (`base_sort_key`/`chem_id`) gains a heading-independent
  canonical ring signature as its leading component in place of arrival index `i`; winding
  sense appended only as a content-identical tiebreak.
- Heading-atom selection for substituted eta rings resolves to `min` by
  `Chem.CanonicalRankAtoms(ringmol)` over the ring's constituent atoms, gated so
  `SYMMETRIC_LIGANDS` continues to win where it currently does.
- Both new computations must derive from the mol/topology already available at the site (no
  behavior-changing re-sanitize); if a canonical rank/SMILES cannot be computed, the code
  falls back to today's behavior (fail-safe, never a silent reorder of unrelated fragments).

### 5.5 Dependencies
- RDKit (`Chem.CanonicalRankAtoms`, `Chem.MolToSmiles`) — already a hard dependency.
- No new libraries.

## 6. Negative Constraints (The "Do NOTs")
- **DO NOT** modify `_determine_winding` or `oin.winding.signed_circulation`, or otherwise
  change how the `>`/`<` character is computed.
- **DO NOT** touch the `SYMMETRIC_LIGANDS` forced-heading path or any non-eta / heterogeneous
  fragment ordering — gate all new logic on "substituted/asymmetric eta group."
- **DO NOT** delete, weaken, un-skip, or auto-substitute
  `test_haptic_face_r2_geometric_fallback_never_auto_substituted`.
- **DO NOT** shift any existing golden: `test_regression_stability` (cisplatin, transplatin,
  cis-PtCl2(en), fac/mer-Ir(ppy)3), the ferrocene / ansa-metallocene haptic goldens, and all
  `tests/candidate_outputs/*` must stay byte-identical.
- **DO NOT** let winding sense enter the RC1 sort key except as the last-resort tiebreak
  between content-identical rings.

## 7. Risks & Mitigation
- **R1 — shared `base_sort_key` blast radius.** `base_sort_key` feeds all fragments, so a
  naive change could reorder unrelated ligands and shift other goldens. → **Mitigation:** gate
  the new identity on eta/haptic groups; keep non-eta identity byte-identical; the full
  existing-golden set is a hard gate (US-003.3).
- **R2 — canonicalization masking a real reflection (CRITICAL).** If winding sense leaked into
  the fragment/heading identity for content-distinct rings, or if the heading rule
  incorporated geometry, a genuinely reflected face could be normalized into a match. →
  **Mitigation:** winding computation frozen (§5.1); winding used in the RC1 key only to
  tiebreak content-identical rings; RC2 heading is topology-only. The R2 skip stays present
  and meaningful as the visible guard (US-004).
- **R3 — canonical rank/SMILES not computable at the site.** The exact mol object at RC1/RC2
  may not be directly rankable without a round-trip. → **Mitigation:** fail-safe fallback to
  current behavior when the signature cannot be computed; verified by the golden gate.

## 8. Success Metrics
- `test_haptic_face_golden_match` un-`expectedFailure`d and PASSES byte-for-byte against
  `tests/candidate_outputs/Ferrocene-halide-face_oin.txt`.
- `uv run python -m unittest discover tests/unit` → OK, **skipped=3, expected failures=0**.
- `uv run python -m unittest discover tests` → OK.
- Every pre-existing golden and `test_regression_stability` string byte-identical.
- The per-ring content-anchored tests (`test_haptic_face_per_ring_flip_inverts_only_that_ring`,
  `test_haptic_face_two_branch_coverage`) still pass; the R2 skip still present and meaningful.
