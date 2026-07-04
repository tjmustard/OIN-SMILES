# MiniPRD: Zone-A P Stereocenter — SPL Dummy-Metal Embed (MiniPRD-C)
**Hypergraph Node ID:** atom_molassembler_adapter
**Parent Node:** mod_generation (`src/oinsmiles/generation/molassembler_adapter.py`)
**Parent SuperPRD:** `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` (v1.0.0)
**Seed decision:** `spec/worklog/SPL-P-enforcement-decision.md` (Option E, `/hyper-consult-cto` 2026-07-03)
**Design brief:** `spec/worklog/SPL-P-enforcement-design-brief.md`
**Execution order:** THIRD in Phase 4 — depends on MiniPRD-A (`_OIN_CIPCode_LP` /
`[P@]` encode contract) and MiniPRD-B (`_verify_zone_a_p` enforcement loop),
both DONE. Fixes the one-sided-SPL correctness bug MiniPRD-B honestly warns on.
**Executor tier:** Sonnet.

## 1. The Confidence Mandate
**Agent Instruction:** Before generating any plans or writing code, analyze this
document and output a Confidence Score (1-10). If below 9, list strictly the
clarifying questions needed to reach 10.

Context an executor MUST load before touching code:
- `generation/molassembler_adapter.py` — `_stitch_fragment` (def ~line 1319; the
  `len(binding_idxs) == 1` monodentate branch at ~1457-1485 with its
  centroid-based orientation heuristic; the `seed` param added by MiniPRD-B), and
  the Zone-A-P verify/re-embed loop (~1916-1982).
- `generation/molassembler_adapter.py` — `_zone_a_p_expected_labels` (~341),
  `_verify_zone_a_p` (~413), and the `normal_frag_meta` records that feed them
  (`frag_smiles`, `binding_idxs`, `atom_start`, `atom_end`).
- `core/chirality.py` — `_build_dummy_metal_copy` (82), `_lp_cip_label` (175),
  `_metal_present_cip_label` (200), and `PseudoAtomStrategy` (22,
  `PSEUDO_ATOMIC_NUM = 0`, wildcard `*`) — the dead-code device this MiniPRD
  revives on the generation side.
- `spec/worklog/SPL-P-enforcement-decision.md` §"Why E", the co-resident-safety
  sketch, and the TET-non-regression argument (both are load-bearing constraints,
  not commentary).

**Root cause being fixed (from the decision note):** generation embeds the
Zone-A-P fragment **trivalent** (metal-free), so `[P@]`/`[P@@]` describe a
3-coordinate centre ETKDG cannot resolve into a metal-facing handedness. On SPL
the metal-present CIP is fixed by which face the metal lands on, and placement
always picks the same face → only one enantiomer is reachable → MiniPRD-B's
re-embed loop can never converge for the other tag and warns "could not be
enforced." The fix makes P 4-coordinate **before** ETKDG so the tag genuinely
controls the embedded geometry, including the metal-facing face.

## 2. Atomic User Stories
* **US-C1:** As the generation pipeline, for a fragment carrying a Zone-A-P tag,
  I attach a Z=0 wildcard as the P's 4th neighbour **before** ETKDG embeds the
  fragment, so `[P@]` and `[P@@]` embed as true 4-coordinate mirror images and
  the dummy atom's embedded position is a faithful metal-facing reference. The
  dummy is stripped after placement — it must never reach `combined_mol` or the
  written XYZ block.
* **US-C2:** As the placement step, for these fragments I orient the fragment by
  aligning the embedded P→dummy vector onto the real slot direction (replacing
  the centroid heuristic for Zone-A-P fragments only), so the metal lands on the
  tag-determined face. Both enantiomers become reachable.
* **US-C3:** As the enforcement contract (MiniPRD-B), `_verify_zone_a_p` now
  passes on the FIRST embed for both `[P@]` and `[P@@]` on `[Pt_SPL]`; the
  re-embed loop stays as a safety net and, when entered, converges (new seed
  still honours the 4-coordinate tag) rather than exhausting 3 attempts and
  warning.
* **US-C4:** As the co-resident-safety invariant, correcting the target P never
  changes any other stereocentre in the same fragment — guaranteed structurally
  (single ETKDG embed satisfies all tags; only the target P gains a neighbour;
  no reflection, no bookkeeping), not by test alone.

## 3. Implementation Plan (Task List)
- [ ] **Task 1 — Zone-A-P detection at the fragment level.** In `_stitch_fragment`,
      determine whether the fragment carries a Zone-A-P tag by reusing
      `_zone_a_p_expected_labels(frag_smiles)` (non-empty ⇒ this fragment has a
      metal-binding P stereocentre). Gate ALL new behaviour on this being
      non-empty. Tag-free fragments take the existing code path **verbatim**
      (inertness guarantee — Task 8 proves byte-identity).
- [ ] **Task 2 — Attach the dummy before ETKDG.** For a Zone-A-P fragment, before
      the embed, build the fragment mol with a Z=0 wildcard bonded to the
      stereogenic P (SINGLE bond, `NoImplicit(True)` on both dummy and P — the
      same valence-bookkeeping fix `_build_dummy_metal_copy` documents at
      `chirality.py:105-108`, so RDKit perceives 4 real substituents and honours
      the `[P@]` tetrahedral tag through ETKDG). Prefer reviving/reusing
      `PseudoAtomStrategy` (`core/chirality.py:22`) for the attach so encode and
      generate share one wildcard-dummy device; if its current shape doesn't fit
      the embed-time need, factor a small shared helper rather than inlining a
      second copy. Track the dummy's atom index and the P's atom index across the
      embed.
- [ ] **Task 3 — Orient by P→dummy, then strip.** Replace the centroid-based
      monodentate orientation heuristic (`_stitch_fragment` ~1466-1483) *for
      Zone-A-P fragments only* with: translate P onto its target slot position,
      then rotate the fragment about P so the embedded P→dummy vector aligns onto
      the slot direction (`slot_units[0]`). This is the same rigid-orientation
      operation as today, but with an exact metal-facing reference instead of the
      away-from-substituents approximation. Then **delete the dummy atom** from
      both the position array and the returned fragment mol, and re-derive
      `symbols`/indices so downstream `atom_start`/`atom_end`, `binding_idxs`, and
      the returned `(positions, symbols, mol)` are exactly as if the dummy never
      existed. The metal itself is added by the assembly step as today — the dummy
      is only an embed-time scaffold.
- [ ] **Task 4 — Index-integrity audit.** Verify the dummy's presence during embed
      does not shift `binding_idxs` or the `local_p_idx` used by
      `_zone_a_p_expected_labels` / `_verify_zone_a_p` (they key off `frag_smiles`
      and `atom_start`). Attach the dummy such that real-atom indices are stable
      (append the dummy as the highest index, or strip-and-renumber before
      return). Add an assertion/comment pinning this.
- [ ] **Task 5 — Confirm the re-embed loop degrades to a no-op safety net.** With
      Tasks 1-4 correct, `_verify_zone_a_p` returns no mismatches on the first
      embed for both tags on SPL, so the loop (~1936) is not entered. Keep the
      loop unchanged (it still protects against rare ETKDG failures); just verify
      it no longer fires for the acceptance fixtures. Do NOT delete it.
- [ ] **Task 6 — Build SPL fixtures** (see §5). `_MONO_P_SPL` and
      `_MONO_P_CORESIDENT_SPL` inline OINs (self-oracle), plus one genuine 3D XYZ
      fixture (`tests/fixtures/PtCl2-PAMP.xyz`-style, PAMP-type monodentate
      P-stereogenic phosphine on Pt-SPL, built in Avogadro, **not** derived from
      any oinsmiles output — dual-copy to `tests/integration/` per convention).
- [ ] **Task 7 — Acceptance + regression tests** (§5 Tests 1-6). Add to
      `tests/unit/test_zone_a_p_genenforce.py` (same node/file as MiniPRD-B).
- [ ] **Task 8 — Inertness proof.** Assert the tag-free goldens (cisplatin,
      transplatin, cis-PtCl2(en), ferrocene, fac/mer-Ir(ppy)3) remain
      byte-identical (`test_regression_stability` + a pre/post XYZ byte-diff on at
      least one). The dummy branch must not be entered for any of them.

## 4. The Negative Space (Constraints)
* **DO NOT** apply a mirror or any improper transform to a fragment — Option R is
  rejected. The correction is a 4-coordinate ETKDG embed; handedness comes from
  the tag, never from reflecting placed coordinates (Resolved B2/B3; decision
  note §"Why not R").
* **DO NOT** change the embed or orientation path for tag-free fragments — gate
  strictly on `_zone_a_p_expected_labels(frag_smiles)` non-empty. Every existing
  golden must stay byte-identical (decision note, TET-non-regression guard #1).
* **DO NOT** let the Z=0 dummy leak past `_stitch_fragment` — it must be stripped
  from positions AND the returned mol before assembly; no stray `*`/Z=0 atom in
  `combined_mol` or the XYZ block. (Red-team's #1 concern.)
* **DO NOT** reimplement the dummy-metal / CIP machinery — reuse
  `PseudoAtomStrategy` / `_build_dummy_metal_copy`'s wildcard-dummy + valence
  fix and the `_lp_cip_label` / `_metal_present_cip_label` recipes from
  `core/chirality.py`. Encode and generate share one convention.
* **DO NOT** regress TET — `[P@]`→R, `[P@@]`→S on the tetrahedral `_MONO_P`
  fixtures must still hold; TET and SPL take the identical new path (decision
  note, guard #2).
* **DO NOT** assert on raw 3D coordinates — assert on derived metal-present CIP
  (`_metal_present_cip_label`), the deterministic oracle over non-deterministic
  ETKDG.
* **DO NOT** touch the DG-fallback path — E is template-path only; DG still warns
  honestly via `_warn_zone_a_p_fallback` (out of scope, decision note §DG).
* **DO NOT** widen scope to Zone-A **N** (RDKit clears trivalent `[N@]`; deferred)
  or to bidentate incompatible-bite routing (TASK-31, separate).

## 5. Integration Tests & Verification
* **Test 1 (Deterministic — SPL enantiomer discrimination) [ACCEPTANCE]:**
  generate from `[Pt_SPL]` `[P@]` and from `[Pt_SPL]` `[P@@]` (same ligand
  skeleton, e.g. `c1ccccc1[P?]{0}(CC)C`, remaining SPL slots filled) →
  `_metal_present_cip_label` gives **opposite** labels, **each matching**
  `_zone_a_p_expected_labels` for its tag, and **neither** emits a
  "could not be enforced" `OINStereoWarning` (run under `-W error::OINStereoWarning`).
  This is the brief's acceptance criterion.
* **Test 2 (Deterministic — TET non-regression):** the existing `_MONO_P_OIN`
  (`[Ni_TET]`) asserts still pass (`[P@]`→R, `[P@@]`→S), zero warnings.
* **Test 3 (Deterministic — co-resident safety on SPL):** `_MONO_P_CORESIDENT_SPL`
  (`[Pt_SPL]`, P stereocentre + directly-bonded carbon stereocentre) → the P's
  metal-present CIP is corrected to the tag AND the co-resident carbon's CIP is
  unchanged vs. the unforced baseline. Extends
  `test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident` to SPL.
* **Test 4 (Deterministic — byte-stable SPL round trip):** XYZ→OIN→XYZ→OIN on the
  new PAMP-type `PtCl2-PAMP` fixture → second OIN byte-identical to the first
  (pinned RDKit); regenerated metal-present CIP matches the original.
* **Test 5 (Deterministic — inertness):** tag-free goldens (cisplatin, transplatin,
  cis-PtCl2(en), ferrocene, fac/mer-Ir(ppy)3) byte-identical; dummy branch not
  entered (assert via a spy/count or a pre/post XYZ byte-diff).
* **Test 6 (Deterministic — loop-as-safety-net):** for Test-1's fixtures,
  `_verify_zone_a_p` returns empty on the first embed (assert the re-embed loop's
  attempt counter stays 0), confirming the fix moved correctness to embed time.
* **Novel:** none — every output here is a deterministic CIP oracle or a
  byte-diff. The one genuine 3D XYZ fixture (`PtCl2-PAMP`) is a **Candidate
  Artifact** requiring human-in-the-loop verification (chemically reasonable
  PAMP-Pt geometry, correct absolute configuration) before promotion to
  `tests/fixtures/`, per the DIPAMP-fixture precedent.

## 6. Acceptance (single-sentence gate)
`[Pt_SPL]` `[P@]` and `[Pt_SPL]` `[P@@]` generate **opposite, correct**
metal-present CIPs with **no "could not be enforced" warning**, plus a
byte-stable SPL round trip — with all existing goldens byte-identical and the
co-resident carbon provably undisturbed.

## 7. Next steps
Optional `/hyper-redteam` → `/hyper-resolve` (this is a decided, tightly-scoped
bug fix seeded by a full consult, so the chain may be short). Then
`/hyper-execute` on Sonnet. After execution: `/hyper-audit`, update
`spec/worklog/NOTES.md` (flip the CONFIRMED-BUG row to DONE) and reconcile
`spec/compiled/architecture.yml`.
