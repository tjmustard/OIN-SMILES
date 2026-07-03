# Process Document: Stereo Phase 3 — Haptic Face Control (`/hyper-resolve`)

**Generated:** 2026-07-03T11:54:42
**Session Focus:** Mediate the Red Team findings for Stereo Phase 3 (Haptic Face Control) and compile the final SuperPRD + MiniPRDs.

## Problem Statement

The OIN→XYZ generation pipeline silently ignores eta-ring winding (`{n>}`/`{n<}`): flipping a
ring's input winding produces byte-identical 3D output, so the prochiral face that coordinates the
metal is uncontrollable from the input OIN. Phase 1 already plumbed the winding signal as far as
`ParsedOIN`/`OINVector`; Phase 3's job is to *consume* it during placement. A Draft PRD and an
adversarial Red Team Report existed; this session ran `/hyper-resolve` to turn those into a
finalized, executable specification.

## Starting State

- Repo at commit `7d90376` (main), clean apart from in-flight worklog edits and untracked
  Phase 3/4 artifacts.
- `spec/active/` held **two parallel sessions' files**: Phase 3
  (`Draft_PRD_StereoPhase3_HapticFace.md`, `RedTeam_Report_StereoPhase3_HapticFace.md`) and the
  parallel Phase 4 Zone-A P-stereocenter work (`Draft_PRD.md`, `RedTeam_Report_ZoneA_P_Encoding.md`).
- The Phase 3 Draft PRD claimed **confidence 9/10**; the Red Team downgraded it to **7/10**,
  raising 8 ranked findings — chiefly an underspecified `signed_circulation` contract, a fixture
  that couldn't exercise the `det +1` justification, and a not-yet-existing golden that risked
  self-blessing a wrong artifact.
- No compiled Phase 3 SuperPRD or MiniPRDs existed yet.

## Approach & Methodology

Spec-driven, following the HACF `/hyper-resolve` state machine: triage highest-severity
collisions first (Phase 1), then group NFR/edge-case gaps under proposed defaults (Phase 2),
confirm Candidate-Artifact routing (Phase 3), then compile and archive (Phase 4). Decisions were
forced as binary/multiple-choice trade-offs via `AskUserQuestion` (max two per turn) rather than
open-ended questions, keeping strictly to vulnerabilities the Red Team raised.

## Steps Taken

1. Read the `/hyper-resolve` SKILL.md, the Draft PRD, and the Red Team Report. Established the
   ranked finding list and confirmed the Red Team's grounding notes (single production callers for
   `_stitch_eta_fragment` and `_determine_winding`; optimiser is rotation-only; golden not yet on
   disk).
2. **Phase 1, batch 1** — asked two forced trade-offs on the two highest-severity findings: the
   `signed_circulation` centroid/ordering/axis contract, and the `det +1` verifiability gap. User
   chose **full closure** for both (contract pinned + behavioral parity test; det assertion +
   chirality-witness fixture).
3. **Phase 1, batch 2** — asked the remaining two Phase-1 findings: seed-masked branch coverage
   and golden origin. User chose **instrumented two-branch coverage** and **pin the golden to the
   hand-verified §5.1 string**.
4. **Phase 2** — grouped findings #5–#8 plus degenerate-input and version-pin gaps into two NFR
   groups (numerical robustness; observability/baselines/pinning) with proposed standard defaults.
   User **approved all** in both groups.
5. **Phase 3** — noted the one non-deterministic output (`Ferrocene-halide-face_oin.txt`) was
   already routed by the golden-pinning decision; the fixture is treated as a Candidate Artifact.
6. **Phase 4** — read the (deprecated) MiniPRD template and the Phase-1 SuperPRD for format
   precedent, then wrote the compiled SuperPRD and two child MiniPRDs.
7. Inspected `archive_specs.py` and discovered it flushes **all** of `spec/active/`. Because the
   directory also held the parallel Phase-4 session's live files, running the script would have
   clobbered them — violating the PRD's own negative constraint. Deviated: manually archived only
   the two Phase-3 files into the script's standard `spec/archive/<timestamp>_<Feature>/` layout,
   leaving Phase 4 untouched.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| `signed_circulation` centers coords internally, pins SMILES order, hosts `n<3` default, fixes axis to metal→centroid outward | Contract-only docstring; document as residual risk | Only full closure prevents the un-centered/mis-ordered "shared" helper from silently recreating the sign-drift bug while passing the §8 grep |
| Add `det(R)≈+1` assertion **and** a chirality-witness ring | det assertion only; defer witness to Phase 4 | The chosen fixture has no stereocenter, so a reflection bug would ship green and detonate in Phase 4; CIP invariance is in-scope for *this operator's* correctness |
| Instrumented two-branch coverage | Flip-only marker test | Seed-42 is deterministic; a marker-flip alone may only exercise one rotation direction, hiding a skip-path sign bug |
| Pin golden to the hand-verified §5.1 string | Machine-emit then human-review | Prevents enshrining a sign bug as "faithful to the wrong answer" |
| Manually archive only Phase-3 files | Run `archive_specs.py` as the skill instructs | The script flushes all of `spec/active/`, which would destroy the parallel Phase-4 session's live work |
| Confidence restored to 9/10 | Leave at Red Team's 7/10 | The three load-bearing findings that justified 7/10 are now fully resolved |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Phase 3 SuperPRD | `spec/compiled/SuperPRD_StereoPhase3_HapticFace.md` | created |
| Helper MiniPRD | `spec/compiled/MiniPRD_SignedCirculationHelper_Phase3.md` | created |
| Correction MiniPRD | `spec/compiled/MiniPRD_HapticFaceCorrection_Phase3.md` | created |
| Phase 3 Draft PRD | `spec/archive/20260703_115343_StereoPhase3_HapticFace/Draft_PRD_StereoPhase3_HapticFace.md` | archived (moved from spec/active) |
| Phase 3 Red Team Report | `spec/archive/20260703_115343_StereoPhase3_HapticFace/RedTeam_Report_StereoPhase3_HapticFace.md` | archived (moved from spec/active) |
| This process document | `spec/process/process_20260703_115442_stereo-phase3-haptic-face-resolve.md` | created |

## Results & Outcomes

All 8 ranked Red Team findings plus the edge-case/NFR gaps now carry documented decisions, every
one taking the full-closure path. A compiled, executable Phase-3 specification exists: one
SuperPRD (system of record, confidence 9/10) and two child MiniPRDs with an explicit dependency
order (helper lands first, correction depends on it). The parallel Phase-4 session's active files
were preserved intact. `spec/active/` now contains only the Phase-4 files.

## How to Reproduce

Prerequisite: repo on `main` at/after `7d90376`; the Phase-3 Draft PRD and Red Team Report present
in `spec/active/` (or restored from `spec/archive/20260703_115343_StereoPhase3_HapticFace/`).

1. Invoke `/hyper-resolve against spec/active/Draft_PRD_StereoPhase3_HapticFace.md and spec/active/RedTeam_Report_StereoPhase3_HapticFace.md`.
2. Answer the Phase-1 forced trade-offs: full closure on the helper contract; det assert + chirality witness.
3. Answer the Phase-1 batch-2 trade-offs: instrumented two-branch coverage; pin golden to the hand-verified string.
4. Approve both Phase-2 NFR groups (numerical robustness; observability/baselines/version-pin).
5. The agent writes the SuperPRD + two MiniPRDs to `spec/compiled/`.
6. **Gotcha (order-dependency):** do NOT run `.agents/scripts/archive_specs.py` while a parallel
   session's files share `spec/active/` — it flushes the whole directory. Archive only the target
   files manually into `spec/archive/<timestamp>_<Feature>/`.

Expected end state: three new files under `spec/compiled/`, two Phase-3 files under the archive
folder, Phase-4 files still in `spec/active/`.

## Patterns & Lessons

- **`archive_specs.py` is directory-global, not selective.** When two `/hyper-*` sessions run in
  parallel and share `spec/active/`, the standard archival step is unsafe — archive by hand and
  keep the script's naming convention for consistency.
- **"Shared helper" is not the same as "no drift."** A single helper that is called with different
  centering/ordering/axis conventions at each site reintroduces the exact bug it was meant to kill,
  and a textual grep for duplicated math will not catch it. The contract (centering, order, default
  guard, axis orientation) must live inside the helper, backed by a fixture-independent behavioral
  parity test.
- **A fixture must be able to exercise the property that justifies the design.** A stereocenter-free
  fixture cannot distinguish a proper rotation (det +1) from a reflection (det −1) — add a witness
  that can, plus a cheap runtime `det` assertion.
- **Deterministic seeds can mask branch coverage.** A fixed embedding seed makes one correction
  branch always-taken; instrument whether the correction fired rather than trusting a marker flip.
- **Goldens should be reasoned, not emitted-then-blessed** — pin to a hand-verified expected string
  so a sign bug fails the test instead of being enshrined.
