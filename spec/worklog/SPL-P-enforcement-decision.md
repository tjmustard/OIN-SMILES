# Decision — fix one-sided Zone-A P enforcement on square-planar

Status: DECISION (output of `/hyper-consult-cto`, Fable, 2026-07-03). Seeds a
MiniPRD via `/hyper-architect`. Supersedes the "TODO — needs a design decision"
row in `NOTES.md`.

Input brief: `spec/worklog/SPL-P-enforcement-design-brief.md`.
Bug of record: `NOTES.md` "CONFIRMED correctness bug" block (commit d44c8d6).

---

## Decision

**Adopt Option E — embed the Zone-A-P fragment with a dummy metal (Z=0) so ETKDG
pins the metal-facing handedness at embed time, symmetric with the encode side's
`_build_dummy_metal_copy`.** Reject P (face-aware placement) and R (reflection
with co-resident protection).

Confidence: 8/10. The mechanism is a direct mirror of code that already works on
the encode side; the residual risk is entirely in the placement/strip plumbing
(mapping the dummy's embedded position onto the real metal slot without
disturbing the current byte-identical output for tag-free fragments), not in the
chirality logic.

---

## Why E — root-cause alignment

The brief's own root cause is the whole argument: on SPL the metal-present CIP is
fixed by which face of the P pyramid the metal lands on, and generation embeds
the ligand **trivalent** (metal-free), so `[P@]`/`[P@@]` describe a 3-coordinate
centre that ETKDG cannot resolve into a metal-facing handedness. Re-seeding
varies the conformer, never the face — hence one reachable enantiomer, hence the
honest-but-wrong "could not be enforced" warning on half the inputs.

E removes the ambiguity at its source: attach a Z=0 wildcard as P's 4th neighbour
**before** ETKDG (SMILES sense: `c1ccccc1[P@]([*])(CC)C`). P is now a genuine
4-coordinate tetrahedral stereocentre; ETKDG honours the tag *including where the
dummy sits*, so `[P@]` and `[P@@]` embed as true mirror images and the dummy's
embedded position is the correct metal-facing reference. Placement then aligns
the P→dummy vector onto the real slot direction and deletes the dummy. Both
enantiomers become reachable, the metal-present CIP is correct on the first
embed, and the re-embed loop reverts to the safety net it already claims to be
(§1936-1968) — a new seed still honours the 4-coordinate tag, so it converges
instead of warning.

This is the exact device the encode side uses (`core/chirality._build_dummy_metal_copy`,
Z=0 on P, valence-bookkeeping fix) and the exact purpose `PseudoAtomStrategy`
(`core/chirality.py:22`, `PSEUDO_ATOMIC_NUM = 0`, wildcard `*`) was built for and
never wired up (D-6). E revives that dead code on the generation side. Encode and
generate end up using one symmetric convention for "represent the metal as a
lowest-priority 4th substituent," which is the property that makes the round trip
trustworthy.

## Why not P (face-aware placement)

P is under-defined for the monodentate case and collapses into E or R. CIP is
**invariant under proper rotation of the whole complex**, so "pick the fragment
orientation whose metal-present CIP matches the target" cannot be satisfied by
rotating a fixed, rigid trivalent embed — every proper reorientation that keeps
the metal on the chemically-mandatory lone-pair face yields the *same* CIP. The
only way to reach the other CIP from a fixed pyramid is to change the pyramid's
handedness, which means either reflecting it (that is R) or re-embedding it with
the tag actually constraining a 4th neighbour (that is E). "Put the metal on the
other face" is not a real alternative: the other face is the substituent-crowded
side, which is both sterically untenable and chemically wrong (the lone pair
binds the metal). P also degrades further for higher denticity, exactly where the
brief flags it. Net: P is not an independent option.

## Why not R (reflection with co-resident protection)

R directly re-introduces the fragility the codebase deliberately avoids. A
whole-fragment reflection inverts *every* stereocentre in the fragment — the
target P plus any co-resident centre (DIPAMP's second P; the co-resident carbon
in the TASK-32 `_MONO_P_CORESIDENT_OIN` fixture). "Un-flip the co-resident ones"
is not a rigid operation: re-inverting a single stereocentre in placed 3D
coordinates means locally reflecting its neighbourhood, which generally forces a
re-embed of that centre anyway — so R ends up doing embed-level work with hand
-rolled bookkeeping on top. It is the one option whose failure mode is silently
corrupting a *different* stereocentre than the one being fixed, and
`test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident`
(`test_zone_a_p_genenforce.py`) exists specifically to keep that from happening.
R optimises for smallest diff at the cost of the safety invariant. Rejected.

---

## Co-resident-safety argument (sketch) — why E is safe by construction

The safety property: correcting the target Zone-A P must never change the
configuration of any other stereocentre that shares its fragment.

Under E the correction happens **entirely inside one ETKDG embed of one
fragment**. The only graph edit is adding a Z=0 wildcard neighbour to the
*target* P; every other stereocentre keeps its own `@`/`@@` tag untouched, and a
single ETKDG call satisfies all tags in the fragment simultaneously (this is the
normal, already-trusted behaviour — ligand `@/@@` passes straight through
embedding, per the TASK-10 finding). There is no reflection and no per-centre
bookkeeping, so there is nothing that *can* invert a co-resident centre. The
dummy is stripped after placement, restoring the exact fragment graph the rest of
the pipeline expects.

This means `test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident`
stays meaningful and passing, and it *extends* naturally to an SPL fixture: build
`_MONO_P_CORESIDENT` on `[Pt_SPL]` (P stereocentre + a directly-bonded carbon
stereocentre) and assert both the P's metal-present CIP is corrected AND the
carbon's CIP is unchanged.

---

## TET-non-regression — how it's guaranteed

TET already reaches both enantiomers (`@`→R, `@@`→S). Two guardrails keep E from
disturbing that:

1. **Scope the dummy embed to fragments that actually carry a Zone-A-P tag.**
   Tag-free fragments (NH3, Cl, Cp, ordinary carbons) keep the existing embed and
   the existing centroid-based monodentate orientation heuristic
   (`_stitch_fragment` §1466-1483) **unchanged** — so every current golden
   (cisplatin, transplatin, ferrocene, fac/mer-Ir(ppy)3, BDPP/BDNN) is
   byte-identical by construction. This is the primary inertness guarantee and it
   should be enforced structurally (branch only entered when
   `_zone_a_p_expected_labels(frag_smiles)` is non-empty), not merely tested.

2. **The dummy embed replaces the centroid heuristic only for Zone-A-P fragments,
   and it is a strictly better metal-facing reference** (an exact embedded
   P→dummy vector vs. an away-from-substituent-centroid approximation). TET and
   SPL take the identical code path; the geometry differs only in the slot
   vectors, which already come from the template. So E fixing SPL cannot leave
   TET behind — they are enforced by the same mechanism.

Regression is proven by the existing `test_zone_a_p_genenforce.py` suite, whose
`_MONO_P_OIN` / `_MONO_P_CORESIDENT_OIN` fixtures are on `[Ni_TET]` today (TASK-32),
plus the new SPL fixtures below. Run both geometries through the same asserts.

---

## Fixture & oracle plan

- **Oracle:** `_metal_present_cip_label` (`core/chirality.py:200`) on the
  regenerated, assembled, DATIVE-bonded 4-coordinate complex — the same recipe
  `_verify_zone_a_p` uses. Expected labels come from the OIN's own `[P@]`/`[P@@]`
  tag via `_zone_a_p_expected_labels` (graph-based `rdCIPLabeler` recompute), so a
  hand-written inline OIN is its own oracle; no OIN-pipeline dependency.
- **New SPL fixtures (build these; independent of the XYZ→OIN pipeline):**
  - `_MONO_P_SPL` — a PAMP-type monodentate P-stereogenic phosphine on Pt/Pd-SPL,
    e.g. `[Pt_SPL].c1ccccc1[P@]{0}(CC)C.[Cl]{1}.[Cl]{2}` (fill the remaining SPL
    slots to a valid 4-coordinate complex). Monodentate on the template path
    (avoids the DIPAMP incompatible-bite DG routing from TASK-30/31), so it
    exercises the enforcing Kabsch path end to end.
  - `_MONO_P_CORESIDENT_SPL` — same, with one P substituent itself a carbon
    stereocentre, for the co-resident assertion.
- **Byte-stable round trip:** a genuine 3D XYZ fixture (PAMP-type, built in
  Avogadro like the DIPAMP fixture, *not* derived from oinsmiles output) →
  `XYZToSMILES` → `OIN3DGenerator.generate()` → `XYZToSMILES` must reproduce
  OIN(1) byte-for-byte. Scope note: byte-stability also depends on non-stereo
  generation fidelity; a *monodentate* phosphine stitches cleanly on the template
  path (unlike bidentate DIPAMP), so this fixture is expected to round-trip
  clean — but the hard stereo gate is the metal-present-CIP asserts, and the
  byte-stable round trip is the fixture-validated goal.

## DG-fallback interaction (scope boundary)

E lives in `_stitch_fragment`'s ETKDG embed on the **template path only** — the
same stage `_verify_zone_a_p` already runs. The DG fallback (`_molassembler_worker`,
and its RDKit-ETKDG sub-fallback) is a separate path that does not run this
enforcement and already warns honestly via `_warn_zone_a_p_fallback` (Phase-4b,
RISK-9; `test_dipamp_dg_fallback_warns_honestly`). **E does not change DG-fallback
behaviour** — a Zone-A-P input that routes to DG still warns "stereo unenforced on
fallback path." Extending the dummy device into the DG path is explicitly out of
scope for this fix; the acceptance fixtures must stay on the template path.

---

## Acceptance test (for the MiniPRD)

1. `[Pt_SPL]` `[P@]` and `[Pt_SPL]` `[P@@]` (same ligand skeleton) generate
   **opposite** metal-present CIPs via `_metal_present_cip_label`, **each correct**
   (matching `_zone_a_p_expected_labels`), and **neither** emits the
   "could not be enforced" `OINStereoWarning` (run under `-W error::OINStereoWarning`).
2. No TET regression: the existing `_MONO_P` / `_MONO_P_CORESIDENT` asserts on the
   tetrahedral fixtures still pass; `test_regression_stability` (cisplatin/
   transplatin/cis-PtCl2(en)/fac+mer-Ir(ppy)3) still byte-identical.
3. Co-resident safety: `_MONO_P_CORESIDENT_SPL` corrects the P CIP while leaving
   the co-resident carbon's CIP unchanged (extends
   `test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident` to SPL).
4. Byte-stable SPL round trip on the PAMP-type XYZ fixture.
5. Inertness: all tag-free goldens (ferrocene, NH3/Cl complexes) byte-identical
   (dummy branch not entered).

---

## Next steps

- `/hyper-architect` → MiniPRD from this note (node
  `atom_molassembler_adapter`; touches `_stitch_fragment` in
  `generation/molassembler_adapter.py` and revives the `PseudoAtomStrategy`
  wildcard-dummy device from `core/chirality.py`). Execute on Sonnet.
- Suggested MiniPRD structure: (1) dummy-attach + embed in `_stitch_fragment`
  gated on Zone-A-P presence; (2) P→dummy → slot alignment replacing the centroid
  heuristic for those fragments; (3) strip dummy pre-assembly; (4) SPL fixtures +
  the 5 acceptance asserts; (5) inertness gate on tag-free fragments.
- Red-team focus: dummy strip completeness (no stray Z=0 atom reaches
  `combined_mol` / the XYZ block); the P→dummy alignment must not perturb the
  Kabsch placement of tag-free fragments; interaction of the dummy atom index with
  `atom_start`/`atom_end` offsets used by `_verify_zone_a_p`.
