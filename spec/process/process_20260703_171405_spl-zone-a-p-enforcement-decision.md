# Process Document: SPL Zone-A P Enforcement — Consult → Decision → MiniPRD

**Generated:** 2026-07-04T00:14:05Z
**Session Focus:** Decide how to fix the one-sided Zone-A P enforcement bug on
square-planar (SPL) transition-metal complexes, and turn that decision into an
executable MiniPRD — using the HACF `/hyper-consult-cto` → `/hyper-architect`
chain.

## Problem Statement

On square-planar complexes, generation of a metal-bound (Zone-A) phosphorus
stereocentre could only ever produce ONE of the two enantiomers: `[Pt_SPL]`
with `[P@]` and `[P@@]` both generated the same metal-present CIP, so ~half of
such inputs emitted the geometrically WRONG enantiomer (honestly warned with
"could not be enforced," but wrong). Tetrahedral (TET) complexes were unaffected.
The fix required a design decision among three competing options, not a quick
patch — that decision, and the resulting MiniPRD, were this session's goal.

## Starting State

- Git HEAD `ad57643` (branch `main`, ahead of origin, clean tree at session
  start). The bug had already been confirmed and documented the same day:
  commit `d44c8d6` recorded it in `spec/worklog/NOTES.md` as a CONFIRMED
  correctness bug, and `ad57643` added `spec/worklog/SPL-P-enforcement-design-brief.md`
  laying out three fix options (E / P / R) plus constraints.
- Stereo Phase 4a (Zone-A P encode) and 4b (gen-enforce) were both DONE: the
  encode side uses a Z=0 dummy-metal device (`core/chirality._build_dummy_metal_copy`)
  to make trivalent P 4-coordinate for CIP labelling; the generation side
  (`molassembler_adapter._verify_zone_a_p` + a bounded re-embed loop) verifies
  and re-embeds on mismatch but deliberately never reflects.
- The gap: the re-embed loop varies the ETKDG *conformer*, but on SPL the
  metal-present CIP is fixed by *placement* (which face of the P pyramid the
  metal lands on), so re-embedding can never flip it. `PseudoAtomStrategy`
  (`core/chirality.py:22`, wildcard `*`, Z=0) existed for exactly this case but
  was dead code.

## Approach & Methodology

Spec-driven, using two HACF phases in sequence: `/hyper-consult-cto` (a
reasoning-heavy CTO persona) to *decide* among options E/P/R against the brief's
hard constraints, then `/hyper-architect` to *turn the decision into a MiniPRD*.
The session ran autonomously (no live user), so the consult made the call from
the brief + code rather than interviewing, and the architect generated the
MiniPRD directly because the decision note was already exhaustive (codebase-first
rule: nothing left to ask). Grounding the decision required reading the actual
enforcement loop and the encode-side dummy-metal machinery before reasoning, so
the trade-off analysis was against real code, not the brief's summary alone.

## Steps Taken

1. Read `spec/worklog/SPL-P-enforcement-design-brief.md` and the first ~700 lines
   of `spec/worklog/NOTES.md` to load the bug, the three options, and the Phase-4
   history. Reason: the brief explicitly forbids jumping to implementation and
   mandates the consult path.
2. Read the rest of `NOTES.md` (the Phase-4b, TASK-30/31/32 log entries) and the
   key source: the Zone-A-P verify/re-embed loop
   (`molassembler_adapter.py` ~1916-1982), the monodentate placement branch of
   `_stitch_fragment` (~1457-1485), and `core/chirality.py`'s
   `_build_dummy_metal_copy` / `_lp_cip_label` / `_metal_present_cip_label`.
   Reason: options E (dummy embed) and P (face-aware placement) hinge on exactly
   how the monodentate path orients the fragment and how the encode-side dummy
   works.
3. Invoked `/hyper-consult-cto`. Key analytical findings that drove the verdict:
   - **P collapses.** CIP is invariant under proper rotation of the whole
     complex, so "pick the orientation whose metal-present CIP matches" cannot be
     satisfied by rotating a rigid trivalent embed — reaching the other CIP
     requires changing the pyramid's handedness (reflect = R) or re-embedding
     with the tag on a real 4th neighbour (= E). Not an independent option.
   - **R re-introduces the exact fragility the codebase avoids** — whole-fragment
     reflection inverts every co-resident stereocentre, and un-flipping one in
     placed 3D isn't rigid, so it collapses into embed-level work with hand-rolled
     bookkeeping. The existing safety test guards against precisely this.
   - **E fixes the root cause the brief names** (trivalent P can't pin the
     metal-facing handedness → make it 4-coordinate), is symmetric with the
     encode side, revives the dead `PseudoAtomStrategy`, and is co-resident-safe
     by construction.
4. Wrote the decision note `spec/worklog/SPL-P-enforcement-decision.md` (verdict E,
   confidence 8/10, with the co-resident-safety sketch, the two-part
   TET-non-regression argument, the fixture/oracle plan, and the DG-fallback
   scope boundary).
5. Updated `spec/worklog/NOTES.md` — flipped the CONFIRMED-BUG status row to
   DECIDED (Option E) and rewrote the "NEXT: consult" block in the bug section to
   record the decision and rationale.
6. Invoked `/hyper-architect`. Read the (deprecated) MiniPRD template and the
   existing `MiniPRD_ZoneA_P_GenEnforce_AUDITED.md` to match house style, node
   ID, and parent SuperPRD. Because the decision note was self-contained, skipped
   the interview and generated the MiniPRD directly.
7. Wrote `spec/active/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` (MiniPRD-C under
   `SuperPRD_StereoPhase4_ZoneA_P.md`, node `atom_molassembler_adapter`): 8 tasks,
   4 user stories, negative-space constraints, 6 tests, executor tier Sonnet.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| Fix via Option E (Z=0 dummy attached before ETKDG) | P (face-aware placement); R (reflection + un-flip co-residents) | E fixes the root cause at embed time, is symmetric with the encode side, revives dead `PseudoAtomStrategy`, and is co-resident-safe by construction |
| Reject P | — | CIP is proper-rotation-invariant, so P can't reach the other enantiomer without becoming E or R; degrades further for higher denticity |
| Reject R | — | Whole-fragment reflection inverts co-resident stereocentres; un-flipping them isn't rigid and re-introduces the exact fragility the safety test guards |
| Gate the dummy branch on `_zone_a_p_expected_labels` non-empty | Apply new orientation to all monodentate fragments | Keeps every tag-free golden byte-identical (inertness) and guarantees TET non-regression structurally |
| Generate MiniPRD directly, skip the interview | Full one-question-at-a-time architect interview | Decision note was exhaustive; codebase-first rule says don't ask what's already resolved; session is autonomous |
| Byte-stable round-trip fixture flagged as human-in-the-loop Candidate Artifact | Auto-derive a fixture from the pipeline | A genuine PAMP-type Pt-SPL XYZ must be built independently (Avogadro) to be a valid oracle, per the DIPAMP-fixture precedent |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Decision note | `spec/worklog/SPL-P-enforcement-decision.md` | created |
| MiniPRD-C (Option E) | `spec/active/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` | created |
| Worklog | `spec/worklog/NOTES.md` | updated (status row + bug-block decision record) |
| This process document | `spec/process/process_20260703_171405_spl-zone-a-p-enforcement-decision.md` | created |

No `src/` or test code was modified this session — it is a design/spec session only.

## Results & Outcomes

A committed-quality design decision (Option E) with a complete rationale, and an
executable MiniPRD ready for `/hyper-execute` on Sonnet. The MiniPRD's acceptance
gate is concrete and testable: `[Pt_SPL]` `[P@]` and `[P@@]` generate opposite,
correct metal-present CIPs (`_metal_present_cip_label`) with no "could not be
enforced" warning, plus a byte-stable SPL round trip, with all existing goldens
byte-identical and the co-resident carbon provably undisturbed. The worklog's
CONFIRMED-BUG row now reads DECIDED and points at the decision note.

## How to Reproduce

Prerequisite: HACF workspace at/after commit `ad57643` on `main`; the design
brief `spec/worklog/SPL-P-enforcement-design-brief.md` and the confirmed-bug
block in `spec/worklog/NOTES.md` present.

1. Read the design brief and `NOTES.md` (bug block + Phase-4/TASK-30..32 log
   entries) for context.
2. Read the load-bearing code before deciding: `molassembler_adapter.py`
   `_verify_zone_a_p` + the re-embed loop (~1916-1982) and `_stitch_fragment`
   monodentate branch (~1457-1485); `core/chirality.py` `_build_dummy_metal_copy`
   (82), `_lp_cip_label` (175), `_metal_present_cip_label` (200),
   `PseudoAtomStrategy` (22).
3. Run `/hyper-consult-cto` seeded by the brief. Expect the E-verdict to fall out
   of two observations: P is proper-rotation-invariant (collapses), R breaks
   co-resident safety. Write the decision to
   `spec/worklog/SPL-P-enforcement-decision.md`.
4. Update `NOTES.md` (status row → DECIDED; record rationale in the bug block).
5. Run `/hyper-architect` on the decision note. Because the note is
   self-contained, generate the MiniPRD directly (no interview). Match the style
   of `spec/archive/MiniPRD_ZoneA_P_GenEnforce_AUDITED.md`; parent it under
   `SuperPRD_StereoPhase4_ZoneA_P.md`, node `atom_molassembler_adapter`. Save to
   `spec/active/MiniPRD_ZoneA_P_SPL_DummyEmbed.md`.

Gotcha / order-dependency: do the code reads (step 2) BEFORE the consult — the
decisive argument against P (CIP proper-rotation-invariance in the monodentate
placement path) is only visible from the actual `_stitch_fragment` orientation
code, not from the brief's prose.

## Patterns & Lessons

- **Ground trade-off consults in the real code, not the brief's summary.** The
  brief presented P as a peer option; reading `_stitch_fragment` revealed it
  collapses into E or R. A consult that reasons only from the brief would have
  under-weighted this.
- **"Symmetric with the existing working side" is a strong tie-breaker.** E won
  partly because the encode side already solved the identical trivalent-P problem
  with a Z=0 dummy — reusing that convention (and reviving `PseudoAtomStrategy`)
  keeps encode/generate coherent and shrinks the round-trip trust surface.
- **Structural guarantees beat test-only ones.** Gating the new embed path on a
  Zone-A-P tag being present makes inertness and TET-non-regression true by
  construction, not merely by a green suite — the MiniPRD states this explicitly.
- **An exhaustive decision note lets `/hyper-architect` skip its interview.** When
  the seed doc already resolves inputs, outputs, scope, node, fixtures, and
  acceptance, the codebase-first rule means the architect should generate
  directly rather than re-ask.
- **Autonomous consult caveat:** with no live user, the CTO persona's "ask
  clarifying questions" default becomes "decide from evidence and state
  confidence" — captured here as an 8/10 with the residual risk (placement/strip
  plumbing, not chirality logic) named explicitly.
