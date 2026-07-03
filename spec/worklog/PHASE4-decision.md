# PHASE4 decision — Zone-A P/N stereocenter encoding

Date: 2026-07-03
Status: DECIDED (CTO consult). Next step: `/hyper-architect` to turn this into
a MiniPRD. Companion brief: `PHASE4-design-brief.md`.

## Decision

**Option A — lone-pair `[P@]`/`[P@@]` convention — for phosphorus. Nitrogen is
explicitly out of scope for Phase 4a** (see "Nitrogen finding" below; Zone-A N
will need an Option-C out-of-band marker if/when a fixture and demand exist).
Option B (wildcard-metal pseudo-atom) is **rejected**.

No OIN grammar change and no version bump: `@`/`@@` inside fragment SMILES is
existing v3.x grammar; Zone-A P tags are new *content*, not new syntax.
`PseudoAtomStrategy` stays dead code (candidate for deletion in the MiniPRD).

## Evidence (RDKit spikes, run 2026-07-03, scine/rdkit env via `uv run`)

Test molecule: `C[P@](CC)c1ccccc1` (methyl/ethyl/phenyl phosphine) and its
`@@` twin. Findings:

1. **String stability (the brief's stated risk — resolved).** Trivalent
   `[P@]` parses, gets `_CIPCode` R/S from the graph, survives
   `AssignStereochemistry(cleanIt=True, force=True)`, survives
   `AddHs`/`RemoveHs`/`SanitizeMol`, and held its CIP label through 20
   randomized `RenumberAtoms` → canonical-SMILES → reparse round-trips.
2. **ETKDG honors the trivalent tag geometrically.** Signed pyramid volume of
   P's three substituents: `[P@]` gave +3.31/+2.57/+4.73/+2.92/+3.74 over five
   seeds; `[P@@]` gave −2.99/−3.90/−3.27/−3.28/−3.22. Opposite, consistent
   sign every time — the primary generation path builds the correct pyramid
   with **zero changes to the embedding step**.
3. **3D perception limitation (shapes the design).**
   `AssignStereochemistryFrom3D` / `AssignAtomChiralTagsFromStructure` return
   nothing for *trivalent* P — but work fine once P has a 4th neighbour
   (tested with a `*` dummy: CIP in == CIP-from-3D for both senses). So all
   geometry→tag derivation and all round-trip verification must happen while
   P is 4-coordinate (metal or dummy present) — which is exactly how
   `CIPAssigner` already operates encode-side.
4. **Nitrogen finding (scope-limiting).** RDKit clears trivalent `[N@]` as
   non-stereogenic (amine inversion): `C[N@](CC)c1ccccc1` →
   `tag=CHI_UNSPECIFIED, CIP=None` after clean. A fragment-SMILES tag on a
   trivalent N **cannot survive any RDKit pass**, ours or a downstream
   consumer's. Option A is therefore P-only by construction. This also kills
   half of Option B: a 4th neighbour on neutral N forces valence 4 → sanitize
   failure unless we lie with `[N@+]`.

## The label convention (the correctness crux)

CIP labels are context-dependent: with the metal present the metal is
(usually) the *highest*-priority substituent; in the fragment the lone pair is
the *lowest*. The R/S letter can legitimately differ between the two views of
the same physical geometry. The MiniPRD must therefore pin one convention:

**The OIN string stores the fragment-local (lone-pair) CIP sense.** Key
equivalence that makes it computable encode-side, while 3D is still available:
a dummy atom (Z=0) is also lowest-priority, and it sits exactly where the
metal was — i.e. where the lone pair points. So:

> lone-pair CIP label of the fragment ≡ CIP label of the 4-coordinate P in a
> copy of the full mol where the metal is swapped to Z=0 and keeps only its
> bond to this P.

Encode-side sketch (replaces the Zone-A `total_deg < 4` **clear** in
`ChiralityRecoveryUtility.recover()`, `core/chirality.py:156-158`):

- `CIPAssigner.assign_all()` (full mol, 3D present): for each metal-bound P,
  build the dummy-metal copy above, `AssignAtomChiralTagsFromStructure` +
  `AssignStereochemistry` (spike 3 proves this works 4-coordinate), and store
  the result as `_OIN_CIPCode_LP` on the real atom. **No label → not
  fragment-stereogenic → no property** — this is the spurious-tag gate: BDPP/
  BDNN-style symmetric phosphines (two identical phenyls) stay tag-free,
  preserving current behaviour.
- `recover()` (fragment mol): Zone-A P *with* `_OIN_CIPCode_LP` → keep the
  copied chiral tag, recompute fragment-local CIP from the tag (spike 1 proves
  this works trivalent), flip on mismatch — i.e. the **same verify-and-flip
  pattern already used for the 4-coordinate zones**, just keyed on the LP
  property. Zone-A P without the property, and all Zone-A N, keep today's
  clearing behaviour. Carbon `@/@@` paths untouched.

## Generation-side consumer (round-trip sketch)

1. **Parse:** `parse_inline_string()` is regex-only (resolution C-2) and
   already preserves `@/@@` → nothing to do; the tag reaches
   `_template_generate`'s fragment SMILES for free.
2. **Embed:** ETKDG builds the correct P pyramid from the tag (spike 2). No
   change to the embedding step.
3. **Placement risk + enforcement:** Kabsch/template placement constrains the
   P *position*, not which face of the pyramid the metal approaches. So add a
   post-assembly **verify-and-reflect**: on the assembled complex (P is
   4-coordinate — 3D perception works, spike 3), recompute the lone-pair-
   convention label via the same dummy-metal trick and compare to the input
   tag; on mismatch, mirror the fragment across the plane of P's three
   substituents and re-place (Phase-3-style machinery). This check IS the
   resurrected Phase-2 ETKDG experiment.
4. **Molassembler fallback path:** unknown whether `from_smiles` respects
   trivalent `[P@]`; MiniPRD investigation item — if not, set an atom
   stereopermutator in `_molassembler_worker`
   (`molassembler_adapter.py:1418`).

## Fixture / oracle plan

- **Primary fixture:** `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` (both P genuine
  CIP centers). Acceptance: `XYZToSMILES().convert()` emits a tag on both
  `P{0}` and `P{1}`; update golden
  `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt`.
- **Lossless oracle:** XYZ→OIN→XYZ→OIN is byte-stable, AND RDKit
  CIP-from-3D on the *regenerated full complex* (4-coordinate P — perception
  valid) matches the original complex's P labels (R,R). Per resolution H-1,
  CIP-from-3D remains the oracle — but it must run metal-present.
- **Enantiomer discrimination:** hand-flip `@↔@@` in the DIPAMP OIN; the two
  generated structures must yield opposite metal-present CIP labels
  (TASK-20-style test, now testable).
- **Negative controls (regression gate):** BDPP/BDNN goldens byte-identical
  (no spurious Zone-A tags); carbon-chirality round-trips (TASK-10 tests)
  stay green.

## Rejected options — why

- **B (wildcard `*` in fragment SMILES):** strictly dominated. Perturbs
  fragment atom counts, slot indices, and binding-atom bookkeeping on both
  pipelines; requires strip/re-insert passes; breaks for N (valence). Spike
  evidence shows A delivers everything B would, with none of this.
- **C (out-of-band per-slot marker):** grammar version bump plus dual-side
  read/write for zero marginal benefit *for P*, now that A is empirically
  validated. **Retained as the designated future path for Zone-A N**, where
  in-fragment tags are impossible (finding 4). Do not build it speculatively.

## Risks / rollback

- Copied-tag parity after metal removal is unverified in the real fragment
  builder (neighbour-order dependent) — that's exactly what the
  verify-and-flip step absorbs; the MiniPRD adds a unit test on the raw
  fragment mol before `recover()`.
- Non-RDKit OIN consumers must accept `[P@]` on trivalent P (standard
  SMILES; Daylight-legal). Documented, not gated.
- Rollback = revert the `recover()` branch + goldens; encoding is
  content-level, so old strings remain parseable throughout.
