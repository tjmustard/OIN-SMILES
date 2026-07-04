# Design brief — fix one-sided Zone-A P enforcement on square-planar

Status: DESIGN CONSULT needed (not a coding task). Seed a `/hyper-consult-cto`
session with this; the output is a DECISION that then becomes a MiniPRD via
`/hyper-architect`. Do NOT jump to implementation from this brief.

## The bug (verified 2026-07-03, Fable — see NOTES.md, commit d44c8d6)

On square-planar (`SPL`) complexes, Zone-A P (metal-bound phosphorus)
stereocenter enforcement can only ever produce ONE of the two enantiomers.
Measured with `_metal_present_cip_label` on the generated mol, ligand
`c1ccccc1[P?]{0}(CC)C` on `[Pt_SPL]`:

| target | generated metal-present CIP | warning |
|---|---|---|
| `[P@]`  | **S** (WRONG — asked for R) | "could not be enforced" |
| `[P@@]` | **S** (right, by luck) | none |
| `[Pt_TET]` `[P@]` | R | none |
| `[Pt_TET]` `[P@@]` | S | none |

TET enforcement works (both enantiomers reachable); SPL is stuck at one.
Impact: SPL is the most common TMC geometry (Pt/Pd, cisplatin analogs), so ~half
of Zone-A-P-on-SPL inputs emit the WRONG enantiomer — honestly warned, but wrong
3D. This undercuts the Phase-4a/4b "enforcement works" claim for SPL.

## Root cause

Enforcement (`molassembler_adapter.py` `_verify_zone_a_p` + loop :1936-1968):
on a detected mismatch it RE-EMBEDS the offending fragment with a new ETKDG
seed, up to 3 times. It deliberately does NOT reflect ("never a mirror/improper
transform (B2/B3)") to avoid inverting co-resident stereocenters (e.g. DIPAMP's
two P atoms). But the metal-present CIP on SPL is fixed by PLACEMENT geometry —
which face of the P pyramid the metal ends up on — not by the embed's internal
chirality that a new seed varies. So re-embedding can never flip it. The only
operations that would: a reflection (excluded), or a placement that respects the
target face.

Context — encode side already solves the analogous problem with a dummy metal
(`core/chirality._build_dummy_metal_copy`, Z=0 on the P) to make the trivalent P
4-coordinate and unambiguous. Generation currently embeds the ligand WITHOUT the
metal (trivalent P → ETKDG can't pin the metal-facing handedness).

## Options to weigh (the consult picks/refines)

**E — Embed with a dummy metal (symmetric with the encode side).** Attach a
Z=0 dummy to the P before ETKDG so the P is 4-coordinate and ETKDG honors the
`[P@]` tag including which face the metal occupies; the dummy's embedded
position then defines the slot/placement. Naturally handles co-resident centers
(ETKDG respects all tags). Reuses the encode-side machinery.
- Pro: fixes the cause at embed time; symmetric with encode; co-resident-safe.
- Con: changes the fragment embed path for Zone-A-P fragments; dummy must be
  stripped and its position mapped to the real slot; interacts with Kabsch.

**P — Face-aware placement.** Keep the metal-free embed, but at placement choose
between the two mirror-related orientations the one whose resulting
metal-present CIP matches the target. For monodentate P the single-atom
alignment leaves orientational freedom to exploit; pick the correct face.
- Pro: localized to placement; no embed change; per-fragment so co-resident-safe.
- Con: must compute target-face orientation reliably; less obvious for
  higher-denticity or when the binding atom isn't the stereocenter.

**R — Reflection with co-resident protection.** Add the excluded reflection,
but after reflecting the offending fragment, re-invert any co-resident
stereocenter back to its original chirality (reflect-all, then un-flip the ones
that shouldn't have moved).
- Pro: directly flips the stuck center; small delta to current loop.
- Con: the "un-flip co-resident" bookkeeping is the exact fragility they avoided;
  error-prone with multiple co-resident centers.

## Constraints / evaluation axes

- MUST NOT regress TET (already correct: `@`→R, `@@`→S) or the carbon-`@/@@`
  paths, or re-introduce spurious tags on symmetric phosphines (BDPP/BDNN).
- MUST preserve co-resident-stereocenter safety — the property
  `test_single_atom_mis_embed_is_corrected_without_mirroring_co_resident`
  guards. Any option must keep that test meaningful (and it should extend to a
  SPL fixture once fixed).
- Oracle = metal-present CIP (`_metal_present_cip_label`) on the regenerated
  4-coordinate complex; round-trip must be byte-stable for a genuine SPL
  P-stereocenter fixture (build one — a monodentate P-stereogenic phosphine on
  Pt/Pd-SPL, e.g. PAMP-type; independent of the OIN pipeline).
- State how each option interacts with the DG fallback path (which does not run
  this enforcement and warns) so the fix's scope is clear.

## Deliverable of the consult

A decision note (`spec/worklog/SPL-P-enforcement-decision.md`): chosen option +
rationale + a sketch of the co-resident-safety argument + the fixture/oracle
plan + how TET-non-regression is guaranteed. THEN `/hyper-architect` → MiniPRD →
Sonnet executor. Acceptance test: `[Pt_SPL]` `[P@]` and `[P@@]` generate
OPPOSITE, correct metal-present CIPs, both without the "could not be enforced"
warning; plus a byte-stable SPL round-trip.
