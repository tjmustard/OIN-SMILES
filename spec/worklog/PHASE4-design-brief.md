# Phase 4 design brief — encoding metal-bound (Zone-A) P/N stereocenters

Status: DECIDED 2026-07-03 — see `PHASE4-decision.md` (Option A for P;
N deferred to a future Option-C marker). Next: `/hyper-architect`.
Use this to seed a `/hyper-consult-cto` or `/hyper-co-research` session. The
output should be a DECISION on representation, which then becomes a MiniPRD via
`/hyper-architect`. Do NOT jump to implementation from this brief.

## The problem (verified 2026-07-03, TASK-20)

OIN silently drops the stereochemistry of any P/N atom bonded directly to the
metal ("Zone A"). Confirmed with `Rh-RR-DIPAMP-Cl2.xyz` (both P atoms are
genuine CIP stereocenters): `XYZToSMILES().convert()` →
`[Rh_SPL].Cc1ccccc1P{0}(CCP{1}(c1ccccc1)c1ccccc1C)c1ccccc1.[Cl]{2}.[Cl]{3}`
— no `@`/`@@` on either P.

Mechanism: `ChiralityRecoveryUtility.recover()` (`core/chirality.py:155`)
clears the chiral tag for any P/N with `total_degree < 4`. A metal binder,
with the metal excluded from the fragment SMILES by OIN's construction, always
has ≤3 fragment-local neighbours → always stripped. The strip happens BEFORE
the branch that consumes `_OIN_CIPCode`.

**Key point: the information is discarded, not lost.** `CIPAssigner` computes
`_OIN_CIPCode` from the intact 3D structure (metal present) and attaches it to
the atom; `recover()` throws it away. So a fix has ground truth to work from.

This is a **format/encoding** question (how does the OIN *string* represent a
metal-bound stereocenter), not a generation-side bug. It supersedes Phase 2:
there is no `@`/`@@` in any input OIN for generation to preserve, so the
generation-side ETKDG experiment (old Phase 2) only becomes testable AFTER
this decision produces a tag that survives the fragment boundary.

## Prior art in-repo

- The `ChiralPNStereocenters` work (archived: `spec/archive/
  SuperPRD_ChiralPNStereocenters.md`, `MiniPRD_ChiralEncoding_AUDITED.md`)
  already considered this and PUNTED on Zone A: resolution C-1 kept
  `PseudoAtomStrategy` only as a "fallback," and Zone-A atoms end up stripped.
  Phase 4 revisits that punt deliberately.
- `PseudoAtomStrategy` (`core/chirality.py:22`) exists for exactly this case
  (backfill the metal as a wildcard 4th substituent, `PSEUDO_ATOMIC_NUM = 0`)
  but is never called — dead code today.

## Representation options to weigh (the consult's job is to pick/refine)

**A. Lone-pair SMILES convention — `[P@]`/`[P@@]` with 3 explicit neighbours.**
SMILES allows tetrahedral chirality on 3-coordinate atoms with a lone pair
(phosphines, sulfoxides); RDKit supports it. Fix = don't strip in `recover()`;
set the tag from the stored `_OIN_CIPCode`, treating the vacated metal position
as the phantom 4th.
- Pro: smallest change; no OIN grammar change; RDKit-parseable.
- Con: must prove the chirality sense relative to (3 subs + phantom) maps
  consistently to (3 subs + metal), and round-trips stably under canonical
  SMILES reordering.

**B. Wildcard-metal pseudo-atom — revive `PseudoAtomStrategy`.**
Insert a `*` (atomic num 0) standing in for the metal so P/N is 4-coordinate;
emit `@`/`@@` against the 4 neighbours incl. `*`.
- Pro: unambiguous 4-coordinate CIP; scaffolding already exists.
- Con: pollutes fragment SMILES with `*`; perturbs fragment atom-count, slot,
  and binding-atom bookkeeping on both encode and generate sides.

**C. Out-of-band annotation — a per-slot chirality marker.**
Carry Zone-A stereo outside the fragment SMILES, symmetric with how the metal
isomer (slot ordering) and winding (`{n>}`/`{n<}`) are already handled.
- Pro: fragment SMILES stays clean/RDKit-standard; decouples from fragment
  CIP ambiguity; consistent with existing OIN design philosophy.
- Con: new OIN grammar (another version bump); both sides need read/write; the
  most work.

## Constraints / evaluation axes

- Project promise is LOSSLESS round-trip — the chosen encoding must survive
  XYZ→OIN→XYZ→OIN and re-derive the same CIP (RDKit CIP-from-3D is the oracle).
- Must not regress the current carbon-`@/@@` behaviour (already correct) or the
  Zone-A CLEARING for genuinely non-stereogenic P/N (must not emit spurious
  tags on symmetric phosphines like the BDPP/BDNN backbone-only cases).
- Should state how the generation side (the deferred Phase-2 ETKDG step) will
  CONSUME whichever encoding is chosen — the decision isn't done until the
  round-trip consumer is sketched.
- Format-version impact: A and B are in-fragment (no/low grammar change); C is
  a version bump (would be the next OIN version after v3.7).

## Deliverable of the consult

A short decision note (append to this file or a new `PHASE4-decision.md`):
chosen option + rationale + the round-trip consumer sketch + a fixture/oracle
plan (reuse `Rh-RR-DIPAMP-Cl2.xyz`). THEN `/hyper-architect` turns it into a
MiniPRD.
