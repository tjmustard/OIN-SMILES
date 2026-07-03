# Process Document: Zone-A P Stereocenter Encoding — Decision & Draft PRD

**Generated:** 2026-07-03T11:23:56-07:00
**Session Focus:** Decide how OIN encodes a metal-bound (Zone-A) P/N stereocenter, then draft the SuperPRD/MiniPRD to implement it (Stereo Roadmap Phase 4).

## Problem Statement

OIN silently drops the stereochemistry of any phosphorus/nitrogen atom bonded directly to the metal ("Zone A"). This breaks the project's core lossless round-trip promise for chiral phosphines (e.g. DIPAMP-type asymmetric-catalysis ligands). The session's job was to *decide* the representation — weighing three candidate encodings — and turn that decision into a formal Draft PRD, without jumping to implementation.

## Starting State

- Git HEAD `7d90376` ("Worklog: Phase 4 design brief + parallel Phase 3/4 plan").
- `spec/worklog/PHASE4-design-brief.md` existed, framing the decision but explicitly marked "DESIGN CONSULT needed (not a coding task yet)". It laid out three options: **A** lone-pair `[P@]` convention, **B** wildcard-metal pseudo-atom, **C** out-of-band per-slot annotation.
- Root cause was already diagnosed (TASK-20, 2026-07-03): `ChiralityRecoveryUtility.recover()` (`core/chirality.py:155-158`) unconditionally clears the chiral tag on any P/N with `total_degree < 4`; a metal binder is always exactly 3-coordinate in the fragment, so it is always Zone A and always stripped — *before* the branch that would consume the stored `_OIN_CIPCode`. `CIPAssigner` computes the correct code from the intact 3D structure, so ground truth exists; `recover()` throws it away.
- `PseudoAtomStrategy` (`core/chirality.py:22`) existed as dead code — scaffolding for Option B, never called.
- Fixture `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` (both P are genuine CIP centers) and placeholder golden `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` were in place.
- `spec/active/` held only a Phase-3 haptic-face draft; no Phase-4 spec yet.

## Approach & Methodology

Spec-driven, following the HACF pipeline: `/hyper-consult-cto` (decide the representation) → `/hyper-architect` (extract requirements into a Draft PRD). The distinguishing methodological choice was **empirical, not paper, decision-making**: rather than reasoning about RDKit's behavior from docs, the CTO consult ran three throwaway RDKit spikes to test Option A's stated risks directly. The spike results reshaped the decision (they scoped nitrogen *out* entirely, which no amount of a-priori reasoning had surfaced).

## Steps Taken

1. **Read the design brief** (`PHASE4-design-brief.md`) plus `core/chirality.py` and `ROADMAP-stereo.md` to ground the three options in the actual strip mechanism and the superseded Phase-2 context.
2. **Invoked `/hyper-consult-cto`.** Adopted the CTO persona and, instead of arguing the options abstractly, wrote and ran a chirality spike (`p3_chirality_spike.py`) against trivalent `[P@]`. Findings: the tag survives sanitize/AddHs/RemoveHs and 20 randomized `RenumberAtoms`→canonical-SMILES round-trips with a stable CIP label — but `AssignStereochemistryFrom3D` returns *nothing* for a trivalent P.
3. **Ran a second spike** (`p3_volume_spike.py`) to test whether ETKDG geometrically honors the tag despite the perception gap. Measured the signed pyramid volume of P's three substituents across 5 seeds: `[P@]` consistently positive, `[P@@]` consistently negative. Also confirmed 3D perception *works* once P is 4-coordinate (dummy `*` attached). Conclusion: the encoding is generation-consumable with zero embed changes; all geometry↔tag work must happen while P is 4-coordinate.
4. **Ran a nitrogen check** — RDKit clears trivalent `[N@]` as non-stereogenic (amine inversion). This was the decisive finding: an in-fragment tag on a metal-bound N *cannot survive any RDKit pass*, so Option A is phosphorus-only by construction.
5. **Wrote the decision note** `spec/worklog/PHASE4-decision.md`: **Option A for P, N deferred to a future Option-C marker, Option B rejected.** Included the lone-pair label convention (the correctness crux), the dummy-metal equivalence that makes it computable encode-side, the generation-side verify-and-reflect consumer sketch, and the fixture/oracle plan. Updated the brief's status line to point at it.
6. **Invoked `/hyper-architect`.** Walked the 5-phase state machine one question per phase, each with a recommended default:
   - Phase 1: confirmed **one SuperPRD, two MiniPRDs** (A encode-side, B generation-side) — because the success oracle is only checkable with the round-trip closed.
   - Phase 2: confirmed the blast radius against `architecture.yml` nodes and the deletion of the dead `PseudoAtomStrategy`.
   - Phase 3: no auth surface; only note is `[P@]` is standard Daylight SMILES.
   - Phase 4: routed the DIPAMP golden through the **Candidate Artifact protocol** — HITL R,R sign-off + an `rdCIPLabeler` warning on CIP conflict (per the user's explicit request).
7. **Generated the Draft PRD** at `spec/active/Draft_PRD.md` with reviewer instructions (§9) naming the exact file to check.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| Option A — lone-pair `[P@]`/`[P@@]` for P | B (wildcard `*` pseudo-atom); C (out-of-band per-slot marker) | Smallest change, no grammar/version bump; spikes proved string stability + ETKDG pyramid fidelity. B perturbs fragment bookkeeping and fails for N; C is a version bump with no marginal benefit for P. |
| Nitrogen out of scope | Treat N symmetrically with P | RDKit clears trivalent `[N@]` (amine inversion) — an in-fragment N tag cannot survive any RDKit pass. C reserved for N if a fixture/demand appears. |
| Store fragment-local (lone-pair) CIP, computed via dummy-metal copy | Store metal-present CIP sense | CIP label is context-dependent; the dummy atom (Z=0) is lowest-priority and sits where the metal/lone-pair was, so the 4-coordinate dummy-metal label equals the trivalent fragment label — and is perceivable (3D works 4-coordinate). |
| Verify-and-flip in `recover()` (not blind trust) | Trust the copied tag directly | Absorbs any atom-order-dependent parity divergence between encode-side and fragment CIP; mirrors the existing 4-coordinate-zone pattern. |
| Two MiniPRDs, A lands first | One combined encode+generate PRD | A is independently green (tag emission + string round-trip); B's oracle needs the closed round-trip. |
| DIPAMP golden = Candidate Artifact, HITL sign-off | Trust `rdCIPLabeler` as fully automated oracle | The golden encodes whatever the code emits (circular); human confirms R,R. RDKit warning stays advisory. |
| Delete `PseudoAtomStrategy` | Keep as fallback | Only purpose was rejected Option B; dead code invites confusion. |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Phase 4 decision note | `spec/worklog/PHASE4-decision.md` | created |
| Phase 4 design brief (status line) | `spec/worklog/PHASE4-design-brief.md` | updated (points at decision) |
| Draft PRD | `spec/active/Draft_PRD.md` | created |
| Chirality spike | scratchpad `p3_chirality_spike.py` | created (throwaway) |
| Pyramid-volume spike | scratchpad `p3_volume_spike.py` | created (throwaway) |
| This process doc | `spec/process/process_20260703_112356_zone-a-p-stereocenter-encoding-decision.md` | created |

## Results & Outcomes

- A committed **decision** (`PHASE4-decision.md`): Option A for P, N deferred, B rejected — backed by reproducible RDKit spike evidence, not assertion.
- A complete **Draft PRD** (`spec/active/Draft_PRD.md`) ready for `/hyper-redteam`, with 7 user stories, a blast-radius map keyed to `architecture.yml` nodes, negative constraints, a risk log, and four open questions (label parity, enforcement placement, molassembler fallback, spurious-tag gate) seeded for adversarial attack.
- Clear HITL instructions: the reviewer checks `tests/candidate_outputs/Rh-RR-DIPAMP-Cl2_oin.txt` against `tests/fixtures/Rh-RR-DIPAMP-Cl2.xyz` to confirm both P are (R,R).
- No production code changed — this was a decision + specification session by design.

## How to Reproduce

Prerequisites: repo at/after HEAD `7d90376`, `uv`-managed env with RDKit (`uv run python`), branch `main`.

1. Read `spec/worklog/PHASE4-design-brief.md`, `src/oinsmiles/core/chirality.py`, `spec/worklog/ROADMAP-stereo.md`.
2. Run `/hyper-consult-cto` with the Phase-4 decision task. Validate Option A empirically with RDKit spikes:
   - Trivalent `[P@]`/`[P@@]`: check `_CIPCode` survives `AssignStereochemistry(cleanIt=True)`, `AddHs`/`RemoveHs`/`SanitizeMol`, and repeated `RenumberAtoms`→`MolToSmiles`→`MolFromSmiles`. (Expect: stable.)
   - Embed each enantiomer with `EmbedMolecule` over several seeds; compute signed volume `v = (n1−P)·((n2−P)×(n3−P))`. (Expect: opposite, consistent sign per enantiomer.)
   - Check `AssignStereochemistryFrom3D` on trivalent P (expect: None) vs 4-coordinate P with a `*` dummy (expect: correct CIP).
   - Check trivalent `[N@]` after `AssignStereochemistry(cleanIt=True)` (expect: cleared → N out of scope).
3. Write `spec/worklog/PHASE4-decision.md` capturing option choice, lone-pair convention, dummy-metal equivalence, generation consumer sketch, fixture/oracle plan.
4. Run `/hyper-architect spec/worklog/PHASE4-decision.md`. Answer the 5 phase gates: two-MiniPRD split; blast radius + delete `PseudoAtomStrategy`; no auth; DIPAMP golden as Candidate Artifact with HITL + RDKit warning.
5. Output lands at `spec/active/Draft_PRD.md`.
6. Next: new conversation → `/hyper-redteam`.

Gotchas / order-dependencies:
- Run the spikes *before* committing to Option A — the N-clearing result is only discoverable empirically and it changes scope.
- All CIP geometry↔tag operations must be done with P 4-coordinate (metal or dummy present); trivalent 3D perception silently returns nothing.
- HACF requires the Red Team to run in a *fresh* conversation (context isolation) — do not continue this thread into `/hyper-redteam`.

## Patterns & Lessons

- **Spike before you spec.** Three ~30-line RDKit scripts converted a paper trade-off into an evidence-backed decision and surfaced the nitrogen scope-limit that reasoning alone missed. The confidence-mandate section of the PRD cites the spikes directly.
- **Context-dependent CIP needs a pinned convention.** For metal-bound stereocenters the R/S letter differs between metal-present and fragment views; the "dummy atom sits where the lone pair points" equivalence is the trick that makes the fragment label computable while 3D is still available.
- **Design around perception gaps, don't fight them.** RDKit won't perceive trivalent-P chirality from 3D, so every derive/verify step is arranged to run 4-coordinate — a constraint that also dictated where generation-side enforcement must live.
- **Candidate Artifact protocol prevents circular goldens.** A stereochemistry golden emitted by the code under test can't validate that code; HITL sign-off plus an independent `rdCIPLabeler` oracle breaks the loop.
