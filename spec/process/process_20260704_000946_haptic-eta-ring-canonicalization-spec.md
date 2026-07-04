# Process Document: Haptic Eta-Ring Canonicalization — Spec Compilation (architect → redteam → resolve)

**Generated:** 2026-07-04T07:09:46Z
**Session Focus:** Turn the haptic-canonicalization design brief into a red-teamed,
executable MiniPRD on hypergraph node `atom_oin_aligner` — without writing any `src/` code.

## Problem Statement

The XYZ→OIN encoder produces a canonical string that is *not invariant* across two
physically-equivalent inputs of a substituted-eta-ring complex (a hand-built fixture vs. the
ETKDG-generated structure of the same molecule). This surfaces as the sole remaining
`@unittest.expectedFailure` in the suite, `test_haptic_face_golden_match`
(`Ferrocene-halide-face.xyz`), where a round-trip re-encode differs from the pinned golden by
two *winding-preserving* relabelings (ring fragment order swaps; the winding marker drifts to a
different ring atom). The goal of this session was to convert a pre-written design brief into a
finalized spec chain (SuperPRD + MiniPRD) ready for a Sonnet executor, hardening it against the
one thing that would make the fix dangerous: masking a genuine reflection.

## Starting State

- Git `HEAD` = `014612a` ("Worklog: diagnose haptic-face golden-match as encoder
  canonicalization"). Branch `main`, ahead of origin, not pushed (standing instruction).
- Suite green except the one expected failure: `discover tests/unit` → 119 OK, skipped=3,
  expected failures=1 (`test_haptic_face_golden_match`).
- Pre-existing inputs authored in prior sessions:
  - `spec/worklog/haptic-canonicalization-design-brief.md` — the diagnostic (measured diff,
    two root causes pinned to `src/oinsmiles/utils/oin_aligner.py`, the "must stay a pure
    winding-preserving relabeling" correctness trap, and non-regression constraints).
  - `spec/worklog/NOTES.md` — the multi-session worklog; "HAPTIC CANON" task row = TODO.
  - `tests/candidate_outputs/Ferrocene-halide-face_oin.txt` — the oracle golden.
- `spec/active/` empty; `spec/compiled/` held three prior SuperPRDs (StereoPhase1/3/4) +
  `architecture.yml`.

## Approach & Methodology

Spec-driven, following the HACF Phase-1 chain end-to-end in a single session at the user's
explicit request: `/hyper-architect` → `/hyper-redteam` → `/hyper-resolve`. The governing
principle was **codebase-first**: rather than trust the design brief's file:line claims, read
the actual code at each phase and let the code adjudicate. That discipline is what turned a
routine "write the PRD" task into catching two brief-level errors before any code was written.

## Steps Taken

1. **Read the inputs.** Read the design brief and `spec/worklog/NOTES.md` (paged, ~1300 lines).
   Established the two root causes: RC1 fragment order = arrival order; RC2 heading atom =
   geometry-dependent for substituted eta rings. Noted the hard safety constraint: winding
   character must keep coming from `_determine_winding`/`signed_circulation`.

2. **`/hyper-architect` (codebase-first exploration).** Read the pinned code directly:
   `oin_aligner.py` `_reduce_hapticity` (`:230-317`), the homogeneous sort + heading selection
   in `_permute_and_serialize` (`:440-634`), `_determine_winding` (`:636-679`),
   `SYMMETRIC_LIGANDS` (`:156-166`), the architecture node `atom_oin_aligner`, the two key
   tests (golden match + the R2 skip), and the golden file. Surfaced the two genuine
   human-judgment decisions (the canonical *rules*) via a single `AskUserQuestion`:
   - **RC2 heading rule** → RDKit canonical atom rank (lowest = heading). *(user accepted)*
   - **RC1 order identity** → canonical ring SMILES, winding sense as tiebreak only. *(user accepted)*
   Wrote `spec/active/Draft_PRD.md` (SuperPRD schema).

3. **`/hyper-redteam` (adversarial, code-grounded).** Before writing the report, verified the
   PRD's own claims against the source. Three checks changed the outcome:
   - Grepped `base_sort_key` / `"key"` usage → it is assigned at `:308` and **never read**
     (dead code). The live identity is `chem_id` at `:313`, which only permutes *within* an
     equal-`chem_id` bucket. → **RT-1 (BLOCKING):** the brief and Draft PRD named the wrong
     lever; the serialized order comes from `"rank": i`, and content-distinct rings never
     compare. A global re-rank inverts cisplatin (`N` < `[Cl]`).
   - Read `oin/winding.py::signed_circulation` → it is single-edge
     (`cross(v_star, v_next)·axis`). → **RT-4 (BLOCKING):** the winding *character* can depend
     on the heading atom, so RC2's "pure relabeling" claim was unproven.
   - Confirmed no RDKit mol is in the `lig` dict (only `smiles` + coords) → **RT-2/RT-3**
     (mol/index provenance, unstable `lig["smiles"]`, bond-perception divergence). Plus
     **RT-5** (missing hard tests + an incomplete key for identical-content/same-winding rings).
   Wrote `spec/active/RedTeam_Report.md` with a findings ledger.

4. **`/hyper-resolve` (forced trade-offs).** Presented the two BLOCKING findings via
   `AskUserQuestion`:
   - RT-1 → **scoped eta-only rank swap** (permute only same-mass eta fragments among the rank
     slots they already occupy; non-eta ranks + metal rank-0 untouched). *(user accepted)*
   - RT-4 → **prove-and-assert, winding math FROZEN** (start-invariance made load-bearing via a
     per-star test + a live reflection test). *(user accepted)*
   RT-2/RT-3/RT-5 folded in as engineering directives (no user trade-off needed).

5. **Compiled and archived.** Wrote `SuperPRD_EtaRingCanonicalization.md` (v1.0.0, RT
   disposition log) and `MiniPRD_EtaRingCanonicalization.md` (10 tasks incl. a Phase-0
   precondition probe; 5 verification tests) to `spec/compiled/`. Ran the archival script —
   first attempt failed (`python: command not found`); re-ran with `uv run python
   .agents/scripts/archive_specs.py EtaRingCanonicalization` → flushed `spec/active/` to
   `spec/archive/20260704_000832_EtaRingCanonicalization/`.

6. **Worklog update.** Appended a dated Log entry to `spec/worklog/NOTES.md` and updated the
   "HAPTIC CANON" task row from TODO → SPEC READY with the next step
   (`/hyper-execute MiniPRD_EtaRingCanonicalization.md`).

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| RC2 heading = lowest `Chem.CanonicalRankAtoms(breakTies=True)` | Highest-CIP-priority substituent; canonically-first substituent SMILES | Purely topological, order-invariant, already in-pipeline; smallest delta from the existing `SYMMETRIC_LIGANDS` "lowest local index" rule |
| RC1 order = canonical ring SMILES, winding sense tiebreak only | Canonical SMILES only; sorted substituent multiset + winding | Stronger identity than a multiset; winding confined to content-identical ties keeps the safety property (for the halide fixture the rings differ, so winding never enters the key) |
| RC1 mechanism = scoped eta-only rank swap | Global re-rank with golden-preserving tiebreak | Global re-rank inverts cisplatin (`N` < `[Cl]`); scoped swap protects every non-eta golden by construction |
| RC4/winding = prove-and-assert, keep math FROZEN | Rewrite character to be explicitly heading-independent | Rewrite unfreezes `signed_circulation`'s shared contract with the generation side; a per-star assertion + reflection test proves invariance for this construction without touching the math |
| Chain all three phases in one session | One phase per conversation (CLAUDE.md context-isolation guidance) | User explicitly requested the full chain; honored the explicit instruction |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| SuperPRD | `spec/compiled/SuperPRD_EtaRingCanonicalization.md` | created |
| MiniPRD | `spec/compiled/MiniPRD_EtaRingCanonicalization.md` | created |
| Draft PRD | `spec/active/Draft_PRD.md` → archived | created, then archived |
| RedTeam Report | `spec/active/RedTeam_Report.md` → archived | created, then archived |
| Archive dir | `spec/archive/20260704_000832_EtaRingCanonicalization/` | created (Draft_PRD + RedTeam_Report) |
| Worklog | `spec/worklog/NOTES.md` | updated (Log entry + task row) |

No `src/` code or test code was modified. No git commit was made.

## Results & Outcomes

A complete, executable spec chain for the eta-ring canonicalization fix now exists in
`spec/compiled/`, hardened by an adversarial pass that corrected two errors carried in the
original design brief (the dead `base_sort_key` lever; the unproven winding start-invariance).
The MiniPRD encodes a fail-closed guardrail: a **Phase-0 precondition probe** that halts the
executor if the lowest-canonical-rank atom does not match the golden's marked atom — explicitly
forbidding a "just re-pin the golden" shortcut. Acceptance is unchanged and objective:
un-`expectedFailure` `test_haptic_face_golden_match` → byte-for-byte pass; `discover tests/unit`
→ skipped=3, expected failures=0; every existing golden byte-identical.

## How to Reproduce

Prerequisites: repo on `main` at `014612a` (or later with `spec/compiled/` intact), `uv`
installed, RDKit available via `uv sync`.

1. Read `spec/worklog/haptic-canonicalization-design-brief.md` and `spec/worklog/NOTES.md`.
2. `/hyper-architect` seeded by the brief (node `atom_oin_aligner`). When it asks the two
   design questions, answer: RC2 = RDKit canonical atom rank; RC1 = canonical ring SMILES with
   winding tiebreak. Output: `spec/active/Draft_PRD.md`.
3. `/hyper-redteam`. **Do the code checks yourself** — grep `base_sort_key`/`"key"` usage
   (it is dead), read `oin/winding.py::signed_circulation` (single-edge), confirm no mol is in
   the `lig` dict. Output: `spec/active/RedTeam_Report.md` with RT-1…RT-5.
4. `/hyper-resolve`. Answer the two forced trade-offs: RT-1 = scoped eta-only rank swap;
   RT-4 = prove-and-assert with winding math frozen. Output:
   `spec/compiled/{SuperPRD,MiniPRD}_EtaRingCanonicalization.md`.
5. Archive: `uv run python .agents/scripts/archive_specs.py EtaRingCanonicalization`
   (note: `python` alone is not on PATH — must use `uv run python`).
6. Next (not part of this session): `/hyper-execute
   spec/compiled/MiniPRD_EtaRingCanonicalization.md` (Sonnet), then `/hyper-audit`.

## Patterns & Lessons

- **Codebase-first is what caught the bugs.** Both BLOCKING findings (RT-1, RT-4) came from
  reading the actual source during the red-team pass, not from re-reasoning about the brief. A
  design brief — even a careful, measured one — can pin a fix to dead code (`base_sort_key`) or
  assume an unproven invariant (winding start-invariance). Verify file:line claims against the
  file.
- **Frozen contracts deserve an explicit "don't touch" plus a test that proves the freeze is
  safe.** Rather than rewrite `signed_circulation`, the resolution kept it frozen and made
  start-invariance load-bearing via a per-star assertion — cheaper and it protects the shared
  encoder/generation convention.
- **Encode fail-closed guardrails into the MiniPRD.** The Phase-0 precondition probe (STOP if
  canonical-rank heading ≠ golden's marked atom; do not re-pin the golden) turns a latent
  reviewer worry into an executor-time hard stop.
- **Environment gotcha:** `.agents/scripts/*` must be run with `uv run python`; bare `python`
  is not on PATH in this environment.
