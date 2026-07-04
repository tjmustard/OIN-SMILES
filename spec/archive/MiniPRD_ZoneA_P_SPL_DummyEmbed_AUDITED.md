# MiniPRD: Zone-A P Stereocenter — Dummy-Metal Embed (MiniPRD-C) — COMPILED
**Hypergraph Node ID:** atom_molassembler_adapter
**Parent Node:** mod_generation (`src/oinsmiles/generation/molassembler_adapter.py`)
**Parent SuperPRD:** `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` (v1.1.0)
**Seed decision:** `spec/worklog/SPL-P-enforcement-decision.md` (Option E, `/hyper-consult-cto` 2026-07-03)
**Design brief:** `spec/worklog/SPL-P-enforcement-design-brief.md`
**Red Team input:** `spec/active/RedTeam_Report.md` (2026-07-03) — 7 findings, all resolved below.
**Resolution:** `/hyper-resolve` 2026-07-03, Thomas Mustard — see §0 Resolved Trade-offs Log.
**Execution order:** THIRD in Phase 4 — depends on MiniPRD-A (`_OIN_CIPCode_LP` /
`[P@]` encode contract) and MiniPRD-B (`_verify_zone_a_p` enforcement loop),
both DONE. Fixes the one-sided-SPL correctness bug MiniPRD-B honestly warns on,
**and (per resolution) extends stereo enforcement to all bidentate Zone-A-P,
including incompatible-bite chelates — absorbing TASK-31.**
**Executor tier:** Sonnet.

---

## 0. Resolved Trade-offs Log (`/hyper-resolve`, 2026-07-03, Thomas Mustard)

The Red Team flagged 2 blockers, 2 high, 2 medium, 1 low. Every flag is closed
here; the body below is written to the resolved decisions, not the pre-resolution
draft.

- **C1 — Reuse target for the embed-time attach (Blockers 1 & 2).**
  *Issue:* the draft told the executor to revive `PseudoAtomStrategy`
  (deleted from `src/` — `grep` returns nothing; `chirality.py:22` is now
  `from .constants import TRANSITION_METALS_NUM`) and to reuse
  `_build_dummy_metal_copy` (which *raises* `ValueError("no metal neighbour
  found on P atom")` on a metal-free fragment and *converts* an existing metal
  rather than *attaching* a new dummy). Both pointers are dead for the embed-time
  need.
  *Options:* (A) new dedicated `_attach_dummy_metal(mol, p_idx)` helper;
  (B) generalize `_build_dummy_metal_copy` into one attach-or-convert helper.
  *Resolution:* **A** — a **new** dedicated helper. It is net-new code, NOT a
  revival. It is *symmetric* to `_build_dummy_metal_copy`'s valence-bookkeeping
  fix (SINGLE bond, `NoImplicit(True)` on both atoms) but **adds** a Z=0
  neighbour instead of converting one. The Confidence Mandate below reflects this
  honestly (net-new, not revived).

- **C2 — Denticity scope (Finding 3 + scope expansion).**
  *Issue:* the gate `_zone_a_p_expected_labels(frag_smiles)` non-empty is also
  true for bidentate Zone-A-P fragments (DIPAMP, BDPP), whose branch has no
  dummy-strip → a Z=0 atom would leak into the XYZ.
  *Options:* (A) monodentate-only, route bidentate to DG honest-warning;
  (B) template-placed bidentates only; (C) **all** bidentate incl. incompatible-bite.
  *Resolution:* **C** — extend the dummy-embed to **all** Zone-A-P fragments of
  any denticity, **including incompatible-bite chelates** currently routed to DG
  by commit `ee0b3f0`. This **absorbs TASK-31** and **supersedes** two original
  MiniPRD-C constraints (see §4 OVERRIDE notes) and the parent SuperPRD's
  "template-path only" framing. Rationale: the dummy fixes P *handedness*, not
  *bite fit*; forcing incompatible-bite onto template placement would not solve
  its bite geometry. So the enforcement scaffold must reach the Zone-A-P atom on
  **whichever placement path bite-compatibility selects** — template OR DG — so
  no Zone-A-P case is left silently unenforced. The dummy-strip postcondition
  (C4) closes the leak for every path.

- **C3 — Test-1 oracle convention (Finding 4).**
  *Issue:* the draft asserted `_metal_present_cip_label` (real high-priority
  metal) *equals* `_zone_a_p_expected_labels` (Z=0 lone-pair convention) — two
  legitimately-divergent conventions; the assertion can pass by cancellation or
  falsely fail.
  *Resolution:* **like-for-like LP + opposite-metal-present.** The enforcement
  assertion compares `_lp_cip_label(dummy copy)` against
  `_zone_a_p_expected_labels` — the exact same-convention pair `_verify_zone_a_p`
  already trusts (SuperPRD-B1, "never cross-convention"). Enantiomer
  *discrimination* is asserted separately: `_metal_present_cip_label` gives
  **opposite** labels for `[P@]` vs `[P@@]`. **Never** assert cross-convention
  equality.

- **C4 — NFR / robustness defaults (Findings 5, 6, 7 + Test-5 strength).**
  *Resolution:* standard defaults **approved as a group** —
  (1) a hard **postcondition** at the end of `_stitch_fragment`: no returned atom
  has `GetAtomicNum() == 0`, and `len(positions)` equals the dummy-absent
  heavy+H count; the one-sentence acceptance gate gains "…and no generated XYZ
  contains a Z=0/`*` atom."
  (2) **Pinned operation order**: `SetNoImplicit(P)` → attach dummy (SINGLE,
  `NoImplicit(dummy)`) → `AddHs` → embed → align → **strip-before-return**; plus
  a **pre-placement parity assert** — `_lp_cip_label` on the embedded, pre-strip
  fragment must match the input `[P@]`/`[P@@]` tag; a mismatch is a **hard
  failure**, never a silent loop-mask.
  (3) **Test 7 (loop-with-dummy)** via the existing `_test_flip_chiral_idx` seam:
  force ≥1 mis-embed so the re-embed loop fires, then assert convergence + no
  Z=0 atom + stable `len(positions)`.
  (4) **Inertness Test 5** proves the dummy branch was never *entered* for
  tag-free goldens via a branch-entry counter/spy — not merely an output byte-diff.

- **C5 — PAMP fixture gating (Finding 4 §5/§6, Candidate Artifact).**
  *Resolution:* **hard merge gate immediately.** Test 4 (byte-stable SPL
  round trip on `PtCl2-PAMP`) is a blocking acceptance criterion now, trusting
  the Avogadro-built fixture's absolute configuration as-is.
  *Accepted risk (logged, RISK-C3):* round-trip stability is pinned to an
  unverified reference geometry; if the fixture's absolute configuration is
  wrong, the gate enforces a wrong-but-consistent round trip. Reviewer sign-off
  is still recorded in `spec/worklog/` post-hoc, but does not block merge.

---

## 1. The Confidence Mandate
**Agent Instruction:** Before generating any plans or writing code, analyze this
document and output a Confidence Score (1-10). If below 9, list strictly the
clarifying questions needed to reach 10.

**Post-resolution honest score: 7/10.** The blockers are closed and the
monodentate path is well-understood, but two areas are net-new and carry genuine
uncertainty the executor MUST surface, not paper over:
- The attach helper `_attach_dummy_metal` is **new code**, not a revival (C1).
- **Bidentate + incompatible-bite (C2) is the real risk surface.** Two open
  sub-questions the executor must resolve *in the MiniPRD, with a spike*, before
  claiming done:
  - **Q-C1 (bidentate orientation):** with two metal-binding P atoms each carrying
    a dummy, orientation is no longer "align one P→dummy vector onto one slot
    direction." Define the rigid alignment as a least-squares (Kabsch) fit of the
    **set** of embedded P→dummy vectors onto the corresponding slot directions
    (`slot_units[binding_rank]`), with the chelate embedded as one fragment so a
    single ETKDG solution satisfies both tags simultaneously.
  - **Q-C2 (incompatible-bite path selection):** the dummy fixes handedness, not
    bite angle. Determine at fragment level whether template placement can honor
    the bite; if not, keep the DG path (`ee0b3f0`) for *placement* but **still run
    the dummy-embed + `_verify_zone_a_p` enforcement** on that path. Enforcement
    must not depend on the placement path chosen.

Context an executor MUST load before touching code:
- `generation/molassembler_adapter.py` — `_stitch_fragment` (def ~line 1319; the
  `len(binding_idxs) == 1` monodentate branch ~1457-1485; the bidentate branch
  ~1487-1519; the `seed` param added by MiniPRD-B), and the Zone-A-P
  verify/re-embed loop (~1916-1982; seed `42 + attempts*1009` at ~1946; the
  index guard at ~1958; the XYZ writer at ~1985).
- `generation/molassembler_adapter.py` — `_zone_a_p_expected_labels` (~341),
  `_verify_zone_a_p` (~413), `_warn_zone_a_p_fallback`, the incompatible-bite → DG
  routing added by `ee0b3f0`, and the `normal_frag_meta` records
  (`frag_smiles`, `binding_idxs`, `atom_start`, `atom_end`).
- `core/chirality.py` — `_build_dummy_metal_copy` (82, the **convert** helper —
  read for the valence-fix recipe at ~152-158 to mirror, do NOT call it here),
  `_lp_cip_label` (175), `_metal_present_cip_label` (200). **`PseudoAtomStrategy`
  / `PSEUDO_ATOMIC_NUM` no longer exist — do not grep for them, do not revive
  them.**
- `spec/worklog/SPL-P-enforcement-decision.md` §"Why E", the co-resident-safety
  sketch, and the TET-non-regression argument (load-bearing constraints).

**Root cause being fixed (from the decision note):** generation embeds the
Zone-A-P fragment **trivalent** (metal-free), so `[P@]`/`[P@@]` describe a
3-coordinate centre ETKDG cannot resolve into a metal-facing handedness. On SPL
the metal-present CIP is fixed by which face the metal lands on, and placement
always picks the same face → only one enantiomer is reachable → MiniPRD-B's
re-embed loop can never converge for the other tag and warns "could not be
enforced." The fix makes P 4-coordinate **before** ETKDG (a Z=0 dummy as the 4th
neighbour) so the tag genuinely controls the embedded geometry, including the
metal-facing face.

## 2. Atomic User Stories
* **US-C1:** As the generation pipeline, for a fragment carrying one or more
  Zone-A-P tags, I attach a Z=0 wildcard as **each** stereogenic P's 4th
  neighbour **before** ETKDG embeds the fragment, in the pinned order
  `SetNoImplicit(P)` → attach dummy (SINGLE bond, `NoImplicit(dummy)`) → `AddHs`
  → embed, so `[P@]`/`[P@@]` embed as true 4-coordinate mirror images and each
  dummy's embedded position is a faithful metal-facing reference. Every dummy is
  stripped after placement — it must never reach `combined_mol` or the written
  XYZ block (enforced by the §3 Task 6 postcondition).
* **US-C2:** As the placement step, for these fragments I orient the fragment by
  aligning the embedded P→dummy vector(s) onto the real slot direction(s)
  (monodentate: single-vector rotation about P onto `slot_units[0]`; bidentate:
  Kabsch least-squares fit of the vector set onto the corresponding
  `slot_units[binding_rank]`), replacing the centroid heuristic for Zone-A-P
  fragments only, so the metal lands on the tag-determined face(s). Both
  enantiomers become reachable.
* **US-C3:** As the enforcement contract (MiniPRD-B), `_verify_zone_a_p` now
  passes on the FIRST embed for both `[P@]` and `[P@@]` on `[Pt_SPL]`; the
  re-embed loop stays as a safety net and, when entered, converges (new seed
  still honours the 4-coordinate tag). Enforcement runs on **whichever placement
  path is selected** — template OR the DG path for incompatible-bite chelates —
  so no Zone-A-P case is left with the old honest-but-unenforced warning.
* **US-C4:** As the co-resident-safety invariant, correcting the target P never
  changes any other stereocentre in the same fragment — guaranteed structurally
  (single ETKDG embed satisfies all tags; only target P atoms gain a neighbour;
  no reflection, no bookkeeping), not by test alone.
* **US-C5:** As the pre-placement guard, the embedded (pre-strip) fragment's P
  handedness is checked against the input tag via `_lp_cip_label` before
  orientation; a mismatch raises a hard failure rather than being masked as a
  loop non-convergence.

## 3. Implementation Plan (Task List)
- [ ] **Task 1 — Zone-A-P detection at the fragment level (any denticity).** In
      `_stitch_fragment`, determine whether the fragment carries one or more
      Zone-A-P tags by reusing `_zone_a_p_expected_labels(frag_smiles)`
      (non-empty ⇒ this fragment has ≥1 metal-binding P stereocentre). Gate ALL
      new behaviour on this being non-empty. **Do NOT co-condition on
      `len(binding_idxs) == 1`** — per C2 the dummy-embed applies to monodentate
      AND bidentate. Tag-free fragments take the existing code path **verbatim**
      (inertness — Task 8).
- [ ] **Task 2 — New `_attach_dummy_metal(mol, p_idx)` helper (C1).** Add a
      dedicated helper in `core/chirality.py` that **appends** a Z=0 wildcard
      bonded to the stereogenic P (SINGLE bond, `NoImplicit(True)` on both dummy
      and P — the same valence-bookkeeping fix `_build_dummy_metal_copy`
      documents at `chirality.py:152-158`, so RDKit perceives 4 real substituents
      and honours the `[P@]` tetrahedral tag through ETKDG). Append the dummy as
      the **highest** atom index. This is **new code, not a revival**; do not
      call `_build_dummy_metal_copy` (it requires an existing metal neighbour and
      raises on a metal-free fragment). Return the dummy's atom index. For a
      bidentate fragment, call once per metal-binding Zone-A-P atom and track all
      dummy indices.
- [ ] **Task 3 — Pinned operation order + pre-placement parity guard (C4.2).**
      Apply, in order: `SetNoImplicit(P)` on each stereogenic P →
      `_attach_dummy_metal` for each → `AddHs` → embed. Because the dummies are
      appended before `AddHs`, confirm H atoms still land after the dummies OR
      re-derive indices after `AddHs`; assert every dummy is a Z=0 atom at a
      known index post-`AddHs`. After embed, before orientation, assert
      `_lp_cip_label` on the pre-strip embedded fragment matches each input tag;
      a mismatch is a hard failure (raise), not a silent skip.
- [ ] **Task 4 — Orient by P→dummy vector(s), then strip.** Replace the
      centroid-based orientation heuristic *for Zone-A-P fragments only* with:
      monodentate — translate P onto its target slot, rotate about P so the
      embedded P→dummy vector aligns onto `slot_units[0]`; bidentate — Kabsch
      least-squares fit of the set of P→dummy vectors onto the corresponding
      `slot_units[binding_rank]` (Q-C1). Same rigid-orientation operation as
      today, with exact metal-facing reference(s) instead of the
      away-from-substituents approximation. Then **delete every dummy atom** from
      both the position array and the returned fragment mol, and re-derive
      `symbols`/indices so downstream `atom_start`/`atom_end`, `binding_idxs`, and
      the returned `(positions, symbols, mol)` are exactly as if no dummy ever
      existed.
- [ ] **Task 5 — Incompatible-bite path (C2, Q-C2).** For a bidentate Zone-A-P
      chelate, decide at fragment level whether template placement can honor the
      bite (reuse the `ee0b3f0` bite-compatibility check). If **compatible**:
      template-place with dummies (Tasks 2-4). If **incompatible**: keep the DG
      path for *placement* (do not force a distorted template), but **still run
      the dummy-embed + `_verify_zone_a_p` enforcement** on the DG-produced mol so
      the Zone-A-P tag is enforced. Never leave a Zone-A-P atom on an
      unenforced path.
- [ ] **Task 6 — Index-integrity audit + hard postcondition (C4.1, Finding 6).**
      Verify the dummies' presence during embed does not shift `binding_idxs` or
      the `local_p_idx` used by `_zone_a_p_expected_labels`/`_verify_zone_a_p`.
      Emit a **hard `assert`** (not a comment) at the single return choke point of
      `_stitch_fragment`: **no** returned atom has `GetAtomicNum() == 0`, AND
      `len(returned_positions)` equals the dummy-absent heavy+H count. This
      converts "the dummy must not leak" into an enforced postcondition and
      protects the XYZ writer (~1985) from ever seeing a `*`/Z=0 token.
- [ ] **Task 7 — Confirm the re-embed loop degrades to a safety net (and strip it
      too).** With Tasks 1-6 correct, `_verify_zone_a_p` returns no mismatches on
      the first embed for both tags on SPL, so the loop (~1936) is not entered for
      acceptance fixtures. Keep the loop unchanged for rare ETKDG failures — but
      because the loop calls `_stitch_fragment` (seed `42 + attempts*1009`), the
      **dummy attach+strip runs inside the loop too**; the Task 6 postcondition
      guards that path. Do NOT delete the loop.
- [ ] **Task 8 — Build fixtures** (see §5): `_MONO_P_SPL` and
      `_MONO_P_CORESIDENT_SPL` inline OINs (self-oracle); a monodentate 3D XYZ
      fixture `tests/fixtures/PtCl2-PAMP.xyz` (PAMP-type, built in Avogadro, NOT
      derived from any oinsmiles output — dual-copy to `tests/integration/` per
      convention); reuse `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` (compatible-bite
      bidentate) and `tests/fixtures/PdCl2-RR-BDPP.xyz`; add **one
      incompatible-bite bidentate Zone-A-P** 3D fixture (Candidate Artifact) to
      exercise the Task 5 DG-with-enforcement path.
- [ ] **Task 9 — Acceptance + regression tests** (§5 Tests 1-7). Add to
      `tests/unit/test_zone_a_p_genenforce.py` (same node/file as MiniPRD-B).
- [ ] **Task 10 — Inertness proof (C4.4).** Assert the tag-free goldens
      (cisplatin, transplatin, cis-PtCl2(en), ferrocene, fac/mer-Ir(ppy)3) remain
      byte-identical, AND prove the dummy branch is **never entered** for them via
      a branch-entry counter/spy (not just an output byte-diff).

## 4. The Negative Space (Constraints)
* **DO NOT** apply a mirror or any improper transform to a fragment — Option R is
  rejected. Handedness comes from the tag through a 4-coordinate ETKDG embed,
  never from reflecting placed coordinates (Resolved B2/B3; decision note
  §"Why not R").
* **DO NOT** change the embed or orientation path for tag-free fragments — gate
  strictly on `_zone_a_p_expected_labels(frag_smiles)` non-empty. Every existing
  golden must stay byte-identical (TET-non-regression guard #1).
* **DO NOT** let any Z=0 dummy leak past `_stitch_fragment` — every dummy must be
  stripped from positions AND the returned mol before assembly; the Task 6
  postcondition `assert` enforces this on every path, including the re-embed loop.
* **DO NOT** revive `PseudoAtomStrategy` or call `_build_dummy_metal_copy` for the
  attach — those are the deleted/convert-only paths the Red Team flagged. Use the
  new `_attach_dummy_metal` (C1). You MAY mirror `_build_dummy_metal_copy`'s
  valence-fix *recipe* (SINGLE bond + `NoImplicit`), not its body.
* **DO NOT** compare CIP labels across conventions in Test 1 — enforcement asserts
  `_lp_cip_label` == `_zone_a_p_expected_labels` (like-for-like); discrimination
  asserts `_metal_present_cip_label` gives opposite labels (C3).
* **DO NOT** regress TET — `[P@]`→R, `[P@@]`→S on the tetrahedral `_MONO_P`
  fixtures must still hold; TET and SPL take the identical new path (guard #2).
* **DO NOT** assert on raw 3D coordinates — assert on derived CIP labels
  (deterministic oracle over non-deterministic ETKDG).
* **DO NOT** widen scope to Zone-A **N** (RDKit clears trivalent `[N@]`; deferred).
* **~~DO NOT touch the DG-fallback path~~ — OVERRIDDEN by C2.** The dummy-embed +
  `_verify_zone_a_p` enforcement now **also** runs on the DG path for Zone-A-P
  fragments (incompatible-bite chelates). DG still uses its own placement; only
  stereo enforcement is added. Non-Zone-A-P DG fallback is unchanged.
* **~~DO NOT widen scope to bidentate incompatible-bite routing (TASK-31)~~ —
  OVERRIDDEN by C2.** This MiniPRD absorbs TASK-31: all bidentate Zone-A-P,
  including incompatible-bite, gets stereo enforcement. Flip the TASK-31 row in
  `spec/worklog/NOTES.md` to "absorbed into MiniPRD-C" on completion.

## 5. Integration Tests & Verification
* **Test 1 (Deterministic — SPL enantiomer discrimination) [ACCEPTANCE]:**
  generate from `[Pt_SPL]` `[P@]` and from `[Pt_SPL]` `[P@@]` (same ligand
  skeleton, remaining SPL slots filled). **Enforcement assertion:**
  `_lp_cip_label` on each generated fragment's dummy copy equals
  `_zone_a_p_expected_labels` for its tag (like-for-like, C3). **Discrimination
  assertion:** `_metal_present_cip_label` gives **opposite** labels for the two
  tags. **Warning assertion:** neither emits a "could not be enforced"
  `OINStereoWarning` (run under `-W error::OINStereoWarning`). Never assert
  `_metal_present_cip_label == _zone_a_p_expected_labels`.
* **Test 2 (Deterministic — TET non-regression):** the existing `_MONO_P_OIN`
  (`[Ni_TET]`) asserts still pass (`[P@]`→R, `[P@@]`→S), zero warnings.
* **Test 3 (Deterministic — co-resident safety on SPL):** `_MONO_P_CORESIDENT_SPL`
  (`[Pt_SPL]`, P stereocentre + directly-bonded carbon stereocentre) → the P's
  label corrected to the tag AND the co-resident carbon's CIP unchanged vs the
  unforced baseline. Extends
  `test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident` to SPL.
* **Test 4 (Deterministic — byte-stable SPL round trip) [ACCEPTANCE, HARD GATE — C5]:**
  XYZ→OIN→XYZ→OIN on the new PAMP-type `PtCl2-PAMP` fixture → second OIN
  byte-identical to the first (pinned RDKit); regenerated CIP matches the
  original. **Merge-blocking now** (accepted risk RISK-C3: fixture absolute
  configuration trusted as-is; reviewer sign-off recorded post-hoc in
  `spec/worklog/`).
* **Test 5 (Deterministic — inertness) [C4.4]:** tag-free goldens (cisplatin,
  transplatin, cis-PtCl2(en), ferrocene, fac/mer-Ir(ppy)3) byte-identical; the
  dummy branch is proven **not entered** via a branch-entry counter/spy (not a
  byte-diff alone).
* **Test 6 (Deterministic — loop-as-safety-net):** for Test-1's fixtures,
  `_verify_zone_a_p` returns empty on the first embed (the re-embed loop's attempt
  counter stays 0), confirming correctness moved to embed time.
* **Test 7 (Deterministic — loop-with-dummy) [C4.3]:** use the existing
  `_test_flip_chiral_idx` seam to force one mis-embed on an SPL Zone-A-P fixture
  so the loop fires ≥1 attempt, then assert (a) convergence, (b) no Z=0 atom in
  the final mol, (c) `len(positions)` stable. This is the ONLY test that
  exercises the strip on the re-embed path.
* **Test 8 (Deterministic — bidentate + incompatible-bite enforcement) [C2]:**
  DIPAMP (compatible-bite) and the new incompatible-bite fixture → both P labels
  match their tags via `_lp_cip_label` == `_zone_a_p_expected_labels`; the
  incompatible-bite case is enforced on the DG path (no "could not be enforced"
  warning); no Z=0 atom leaks (Task 6 postcondition holds on the DG path too).
* **Candidate Artifacts:** the `PtCl2-PAMP` monodentate fixture and the
  incompatible-bite bidentate fixture are hand-built 3D structures. Per C5 the
  PAMP fixture hard-gates immediately with post-hoc sign-off; the incompatible-bite
  fixture follows the same protocol.

## 6. Acceptance (single-sentence gate)
`[Pt_SPL]` `[P@]` and `[Pt_SPL]` `[P@@]` generate **like-for-like-correct LP
labels and opposite metal-present CIPs** with **no "could not be enforced"
warning**; a byte-stable SPL round trip on `PtCl2-PAMP`; bidentate Zone-A-P
(including one incompatible-bite chelate) enforced on its placement path; **all
existing goldens byte-identical**, the co-resident carbon provably undisturbed,
and **no generated XYZ contains a Z=0/`*` atom**.

## 7. Risks & Mitigation
- **RISK-C1 (bidentate orientation, Q-C1):** two-vector Kabsch alignment may not
  simultaneously satisfy both bite geometry and both tag faces. → single-fragment
  ETKDG embed satisfies both tags by construction; orientation is a rigid fit,
  and Test 8 asserts both labels. Spike the two-vector fit before wiring it in.
- **RISK-C2 (incompatible-bite placement, Q-C2):** pulling incompatible-bite
  chelates back toward template placement could distort geometry. → C2/Task 5
  keeps DG *placement* for incompatible-bite and only adds enforcement; no forced
  template distortion.
- **RISK-C3 (unverified acceptance fixture, C5 accepted):** Test 4 hard-gates on
  an Avogadro fixture whose absolute configuration is trusted as-is. → reviewer
  sign-off recorded post-hoc in `spec/worklog/`; if wrong, regenerate fixture and
  the golden re-derives.
- **RISK-C4 (embed-time parity flip on attach):** adding a 4th neighbour to a P
  RDKit modelled as 3-explicit-plus-lone-pair could reinterpret parity. → the
  pre-placement `_lp_cip_label` parity assert (Task 3) makes any flip a hard
  failure, not a silent wrong-enantiomer emit.
- **RISK-C5 (wildcard embed failure):** `EmbedMolecule` may return −1 on a
  `P–*` (Z=0) bond with no UFF params. → confirm embed succeeds on the
  dummy-bearing mol in a spike; if it fails, the fragment falls to DG-with-dummy
  (Task 5) rather than silently dropping enforcement.

## 8. Next steps
`/hyper-execute` on Sonnet (this compiled MiniPRD closes all 7 Red Team findings).
After execution: `/hyper-audit`, update `spec/worklog/NOTES.md` (flip the
CONFIRMED-BUG row to DONE and the TASK-31 row to "absorbed into MiniPRD-C"), and
reconcile `spec/compiled/architecture.yml` (`atom_molassembler_adapter`;
new `_attach_dummy_metal` in `atom` node under `core/chirality.py`).
