# Process Document: Zone-A P SPL Dummy-Embed — Red Team Resolution (MiniPRD-C)

**Generated:** 2026-07-04T00:38:28Z
**Session Focus:** `/hyper-resolve` pass over MiniPRD-C (`MiniPRD_ZoneA_P_SPL_DummyEmbed.md`) — mediate the Red Team's 7 findings, force the architectural trade-offs, compile the final MiniPRD + update the parent SuperPRD, and archive the active specs.

## Problem Statement

Phase-4 Zone-A phosphorus stereochemistry has a known one-sided-SPL correctness bug: generation embeds the metal-binding P atom **trivalent** (metal-free) before ETKDG, so `[P@]`/`[P@@]` cannot control the metal-facing handedness on square-planar (SPL) centres — only one enantiomer is reachable, and MiniPRD-B's re-embed loop honestly warns "could not be enforced." MiniPRD-C was drafted to fix this by attaching a Z=0 dummy as the P's 4th neighbour before embedding, then Red-Teamed. This session's job was to resolve the Red Team report before `/hyper-execute`.

## Starting State

- Branch `main`, HEAD at `ad57643` ("Worklog: add SPL Zone-A P enforcement design brief").
- `spec/active/` held the two inputs to resolve: `MiniPRD_ZoneA_P_SPL_DummyEmbed.md` (the draft) and `RedTeam_Report.md` (7 findings: 2 blocker, 2 high, 2 medium, 1 low).
- Parent spec `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` was v1.0.0 with MiniPRD-A (encode) and MiniPRD-B (gen-enforce) both DONE; MiniPRD-C unlisted.
- Seed decision `spec/worklog/SPL-P-enforcement-decision.md` (Option E, from a `/hyper-consult-cto`).
- Key pre-existing behaviour in play: commit `ee0b3f0` routed incompatible-bite chelates to DG (placement fidelity), and the seed note deferred TASK-31 (incompatible-bite routing).

## Approach & Methodology

Spec-driven mediation following the `/hyper-resolve` state machine: verify the highest-severity findings against ground truth first, then drive forced trade-offs via `AskUserQuestion` (max 2 per turn, highest-risk first), then compile. Rather than trust the Red Team's claims, the two BLOCKERs were verified directly against `src/` before mediating — both were confirmed real, which changed the framing of the reuse-target question from a doc typo to a genuine architectural fork.

## Steps Taken

1. Read `.agents/skills/hyper-resolve/SKILL.md`, then the two active inputs (`RedTeam_Report.md`, `MiniPRD_ZoneA_P_SPL_DummyEmbed.md`) to load the full finding set and the draft's task list.
2. **Verified both blockers against source** (Bash `grep` + `sed` on `core/chirality.py`), because a resolution built on a stale finding would misdirect the executor. Confirmed: `PseudoAtomStrategy`/`PSEUDO_ATOMIC_NUM` are absent from `src/` (`chirality.py:22` now imports `TRANSITION_METALS_NUM`), and `_build_dummy_metal_copy` opens with `if metal_idx is None: raise ValueError(...)` then *converts* an existing metal to a dummy — the inverse of the *attach* MiniPRD-C needs on a metal-free fragment.
3. **Phase 1 — two forced trade-offs** via `AskUserQuestion`: (Q1) how to build the embed-time attach given both reuse pointers are dead; (Q2) how to handle bidentate Zone-A-P, which also satisfies the gate but had no dummy-strip. User chose **new dedicated helper** and, unexpectedly, **extend to bidentate now**.
4. **Phase 1 continued — two more** trade-offs: (Q3) how far "extend to bidentate" reaches given `ee0b3f0` already routes incompatible-bite to DG and TASK-31 was deferred; (Q4) the Test-1 oracle convention (Finding 4). User chose **all bidentate incl. incompatible-bite** (absorbing TASK-31) and **like-for-like LP + opposite-metal-present**.
5. Made one engineering call inside the C2 expansion and surfaced it explicitly: the dummy fixes P *handedness*, not *bite fit*, so incompatible-bite chelates should **not** be forced onto template placement. Instead the dummy-embed + `_verify_zone_a_p` enforcement runs on **both** template and DG paths for Zone-A-P — overriding MiniPRD-C's original "template-path only" constraint and superseding `ee0b3f0`'s enforcement-free routing.
6. **Phase 2 — grouped NFR defaults + Candidate-Artifact routing** via `AskUserQuestion`: bundled Findings 5/6/7 + Test-5 strength as one "approve defaults" question (approved), and routed the hand-built PAMP fixture (Finding 4 §5/§6). User chose **hard merge gate immediately** for the PAMP fixture — recorded as accepted risk RISK-C3.
7. **Phase 4 — compilation.** Read the MiniPRD template + parent SuperPRD to conform. Wrote the compiled `spec/compiled/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` with a §0 Resolved Trade-offs Log (C1–C5) and every finding closed to the resolved decision.
8. Updated the parent SuperPRD: bumped v1.0.0 → v1.1.0, registered MiniPRD-C in the child list and execution checklist, and added a §5.3.1 SuperPRD-level trade-offs log. Fixed a section-numbering collision (an inserted `5.5` landed before the existing `5.4`/`5.5 Dependencies`) by renumbering the addition to §5.3.1 and correcting the cross-reference.
9. Ran the centralized archival script. `python` was absent (exit 127); fell back to `uv run python .agents/scripts/archive_specs.py ZoneA_P_SPL_DummyEmbed`, which flushed `spec/active/` to the archive folder.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| New dedicated `_attach_dummy_metal(mol, p_idx)` helper (C1) | Revive `PseudoAtomStrategy`; generalize `_build_dummy_metal_copy` into attach-or-convert | Both reuse targets verified dead/inverse; a net-new helper mirroring only the SINGLE-bond + `NoImplicit` valence recipe is honest and avoids widening the encode helper's blast radius |
| Extend dummy-embed to **all** bidentate incl. incompatible-bite (C2) | Monodentate-only + DG honest-warning; template-placed bidentates only | User choice — absorbs TASK-31 so no Zone-A-P case is left silently unenforced |
| Enforcement runs on both template AND DG paths for Zone-A-P (engineering call inside C2) | Force incompatible-bite onto template placement | The dummy fixes handedness, not bite fit; forcing template would distort geometry. DG keeps placement, gains enforcement |
| Test-1 like-for-like LP oracle + opposite-metal-present (C3) | Metal-present as sole normative oracle | LP-vs-LP matches `_verify_zone_a_p` (SuperPRD-B1 "never cross-convention"); cross-convention equality can pass by cancellation or falsely fail |
| Approve NFR defaults as a group (C4) | Modify some; reject | Each is a Red-Team-recommended best practice (postcondition assert, pinned op order + parity guard, Test 7 loop-with-dummy, spy-based inertness) |
| PAMP fixture = hard merge gate now (C5) | Non-blocking until human validation | User choice; residual risk (unverified absolute config) logged as RISK-C3 with post-hoc sign-off |
| Confidence 7/10, not the seed's 8 | Keep 8/10 | Attach helper is net-new and bidentate/incompatible-bite are net-new (Q-C1, Q-C2 open); honesty over optimism |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Compiled MiniPRD-C | `spec/compiled/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` | created (resolved) |
| Parent SuperPRD | `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` | updated (v1.1.0: child list, §5.3.1 log, checklist) |
| Draft MiniPRD-C | `spec/active/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` | archived |
| Red Team report | `spec/active/RedTeam_Report.md` | archived |
| Archive folder | `spec/archive/20260703_173152_ZoneA_P_SPL_DummyEmbed/` | created (holds the two archived files) |

## Results & Outcomes

All 7 Red Team findings carry a documented decision; `spec/active/` is flushed to `.gitkeep`-only state. The compiled MiniPRD-C is ready for `/hyper-execute` on Sonnet, with two blockers closed against verified ground truth and the scope expansion (absorb TASK-31, supersede `ee0b3f0` for Zone-A-P) recorded in both the MiniPRD and the SuperPRD. Two open sub-questions (Q-C1 bidentate Kabsch orientation, Q-C2 incompatible-bite path selection) are flagged as spike-before-done rather than buried, and the honest confidence is 7/10.

## How to Reproduce

Prerequisite state: branch `main` at `ad57643`; `spec/active/` containing `MiniPRD_ZoneA_P_SPL_DummyEmbed.md` and `RedTeam_Report.md`; `uv` available (the repo uses `uv run python`, plain `python` is not on PATH).

1. `/hyper-resolve` — reads the SKILL, the Red Team report, and the draft MiniPRD.
2. Before mediating, verify the blockers: `grep -rn "PseudoAtomStrategy\|PSEUDO_ATOMIC_NUM" src/` (expect empty) and inspect `_build_dummy_metal_copy` in `core/chirality.py` (expect a `raise ValueError` on missing metal neighbour). This reframes Finding 1/2 from typo to real fork.
3. Answer the four Phase-1 trade-offs (attach helper; bidentate handling; bidentate reach vs `ee0b3f0`/TASK-31; Test-1 oracle) and the two Phase-2 questions (NFR defaults group; PAMP fixture gating).
4. Compile: write `spec/compiled/MiniPRD_ZoneA_P_SPL_DummyEmbed.md` and update `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` (version bump + child registration + trade-offs log). Watch section numbering when inserting into an existing §5.x block.
5. Archive: `uv run python .agents/scripts/archive_specs.py ZoneA_P_SPL_DummyEmbed` (not bare `python` — exit 127). Log the returned absolute archive path.

Gotcha / order-dependency: verify blockers **before** asking the trade-off questions, or the reuse-target question is framed wrongly. The Test-1 oracle decision (C3) must land before the acceptance gate is restated, since the gate inherits its convention.

## Patterns & Lessons

- **Verify Red Team blockers against source before mediating.** Both blockers here were real, and the verification changed a "fix the doc pointer" into a genuine architectural fork (new helper vs generalize). A stale finding accepted at face value would have misdirected the executor.
- **When a user's forced-trade-off choice expands scope, name the second-order engineering consequence rather than silently maximizing.** "Extend to bidentate incl. incompatible-bite" doesn't mean "force incompatible-bite onto template" — the dummy addresses handedness, not bite fit — so enforcement was routed onto both placement paths and the two overridden constraints were called out explicitly.
- **Overrides belong in the negative-space section, struck-through with an OVERRIDE note**, so the executor sees exactly which prior DO-NOTs no longer apply and why (`ee0b3f0`, TASK-31).
- **Section-numbering hygiene:** inserting a heading into an existing numbered block (§5.x) can shadow a later same-numbered section; renumber to a sub-section (§5.3.1) to keep ordering monotonic.
- **Environment note:** this repo has no bare `python`; use `uv run python` for `.agents/scripts/*`.
