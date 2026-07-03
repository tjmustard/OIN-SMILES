# Process Document: Zone-A P Stereocenter Encoding — /hyper-resolve Session

**Generated:** 2026-07-03T11:57:10-07:00
**Session Focus:** Mediating the Red Team report against the Phase 4 (Zone-A P encoding) Draft PRD and compiling the final SuperPRD + two MiniPRDs.

## Problem Statement

The Stereo Roadmap Phase 4 Draft PRD (`Draft_PRD.md`, v0.1.0) — which specifies how OIN will stop silently dropping the stereochemistry of metal-bound phosphorus atoms — had been through `/hyper-redteam`, producing a report with ten consolidated findings (B1–B10), three of them blockers. Per the HACF Phase 1 flow, every finding needed a documented human decision before the spec could be compiled into executable MiniPRDs.

## Starting State

- Branch `main` at commit `7d90376` ("Worklog: Phase 4 design brief + parallel Phase 3/4 plan").
- `spec/active/` held four files: `Draft_PRD.md` + `RedTeam_Report_ZoneA_P_Encoding.md` (this feature) and `Draft_PRD_StereoPhase3_HapticFace.md` + `RedTeam_Report_StereoPhase3_HapticFace.md` (a parallel, then-unresolved Phase 3 draft).
- `spec/compiled/` held `architecture.yml` and `SuperPRD_StereoPhase1_Winding.md` — establishing the per-feature SuperPRD naming convention (the old monolithic `SuperPRD.md` no longer exists).
- The Red Team verdict: core design (Option A, lone-pair `[P@]` convention, verify-and-flip) sound and well-evidenced, but with one internal contradiction (cross-convention warning oracle), one dependency on nonexistent machinery (Phase-3 reflection), and one correctness gap (global fragment mirror corrupts co-resident stereocenters).

## Approach & Methodology

Strict execution of the `/hyper-resolve` state machine (`.agents/skills/hyper-resolve/SKILL.md`): Phase 1 triaged the blockers and high-severity items as forced binary/ternary trade-offs via `AskUserQuestion`, two per turn; Phase 2 grouped the Medium/Low findings into a standard-defaults package for single approve/modify/reject; Phase 3 confirmed the Candidate Artifact (HITL golden) routing; Phase 4 compiled the SuperPRD + MiniPRDs and ran the archival script. Every option was framed as cost vs. risk, with a recommended choice listed first.

## Steps Taken

1. Read the skill, the Draft PRD, and the Red Team report in full to extract the consolidated blocking-items table (B1–B10).
2. **Batch 1 (blockers):** presented B1 (warning oracle compares lone-pair vs metal-present conventions the PRD itself declares incomparable) and B2+B3 as one intertwined decision (mis-embed correction cites nonexistent "Phase-3-style" reflection; a global mirror inverts every co-resident stereocenter). User chose: B1 → dummy-copy oracle (same-convention only); B2/B3 → bounded ETKDG re-embed with a new seed, no mirror ever.
3. **Batch 2 (highs):** presented B4 (dummy-metal copy sanitize can throw on eta-ligand complexes → `convert()` crash regression) and B5 (legacy `AssignStereochemistry` vs `rdCIPLabeler` must be pinned end-to-end). User chose: B4 → degrade to store-nothing + warn, with a CpM(PR₃) regression fixture; B5 → `rdCIPLabeler` authoritative everywhere, guarded + skippable.
4. **Phase 2 (NFR defaults):** tabled B6–B10 as a package — recompute-CIP-per-flip (≤2 passes), exactly-one-metal eligibility guard with the predicate imported from `utils/xyz2mol.py:30` (no duplicate list), 3-attempt re-embed budget with mol=None fallback paths skip+warn, `_LP`-property-before-degree branch precedence + `assign_all()` idempotence, `OINStereoWarning(UserWarning)` with atom index in-message and a `-W error::OINStereoWarning` clean-fixture gate, housekeeping (dangling `architecture.yml` edge, stale `:1418`→`:1489` line ref, constraint reword, `[P@]{0}` parse-adjacency test), byte-stability scoped to the pinned RDKit version. Approved wholesale.
5. **Phase 3 (Candidate Artifact routing):** proposed the amended HITL protocol from Red Team §9 — reviewer told explicitly that the OIN tag encodes the lone-pair sense while the printed rdCIPLabeler table shows the metal-present sense (confirm (R,R) against the metal-present column), execution run emits a per-atom label table + mol-block depiction, sign-off recorded in `spec/worklog/`, duplicate fixture (`tests/integration/` copy) deleted with the golden's provenance line carrying fixture path + SHA. Approved.
6. Read the schema templates (`.agents/schemas/`, both deprecated-but-authoritative), the archive script, and the Phase 1 compiled SuperPRD to match format conventions; noted `spec/compiled/` uses per-feature SuperPRD files.
7. Wrote `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` (v1.0.0, confidence 9/10 with an explicit generation-side deduction, full Resolved Trade-offs Log, four new risks RISK-7..10, amended §9 HITL instructions).
8. Wrote `spec/compiled/MiniPRD_ZoneA_P_Encode.md` (14 tasks) and `spec/compiled/MiniPRD_ZoneA_P_GenEnforce.md` (8 tasks, blocked on A) per the MiniPRD schema.
9. **Archival with a guard:** the archive script flushes ALL of `spec/active/`, which at session start also held the unresolved Phase 3 drafts. Planned a stash-run-restore around the script; on execution discovered the Phase 3 files had already been archived by their own parallel resolve session (`spec/archive/20260703_115343_StereoPhase3_HapticFace/`), so the plain flush was safe. Ran `python3 .agents/scripts/archive_specs.py ZoneA_P_Encoding` (note: `python` is not on PATH; `python3` is).
10. Fixed the stale "SuperPRD (system of record)" pointer in the project memory index to reflect the per-feature convention.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| B1: rdCIPLabeler oracle runs on the same dummy-metal copy that produced `_OIN_CIPCode_LP` | Two separate same-convention checks; dropping the auto-warning | Apples-to-apples comparison; a warning means genuine self-inconsistency, never a legitimate cross-convention difference. Metal-present labels stay print-only for HITL. |
| B2/B3: mis-embed correction = re-embed fragment with new ETKDG seed, ≤3 attempts, no mirror | Build a decision-table mirror utility now; block on Phase 3's reflection landing | Co-resident stereocenters are safe by construction (no improper transform exists to corrupt them); zero dependency on the unresolved Phase 3 draft; ETKDG mis-embeds are rare (spike 2) so this is a safety net, not a hot path. |
| B4: dummy-copy failure degrades to store-nothing + `OINStereoWarning` | Fail loud (propagate the sanitize exception) | An unguarded exception would crash `convert()` for eta-ligand+phosphine complexes that round-trip today — a functional regression worse than the status-quo stereo loss. Gap stays observable via the warning. |
| B5: `rdCIPLabeler` authoritative for every `_OIN_CIPCode_LP` computation/recompute/check | Legacy `Chem.AssignStereochemistry` `_CIPCode` end-to-end | Modern implementation is rules-complete (CIP 4b/5); mixing the two creates a permanent unactionable warning state. Guarded per-call + skippable flag covers the pathological-runtime risk. |
| B6–B10 NFR package approved as tabled | Per-item modification | Red Team's proposed defaults were already conservative and internally consistent. |
| HITL protocol amended (convention statement, label table + depiction, recorded sign-off, fixture dedup) | Table-only; keep original §9 | Without the convention statement a correct string could be rejected (or a wrong one approved) by a reviewer thinking in metal-present CIP. |
| Per-feature SuperPRD file (`SuperPRD_StereoPhase4_ZoneA_P.md`) | Overwriting a monolithic `SuperPRD.md` | Matches the existing `SuperPRD_StereoPhase1_Winding.md` convention in `spec/compiled/`; the monolith no longer exists. |
| Stash-guard around the archive flush (ultimately unneeded) | Run the flush blind; ask the user | The script moves everything in `spec/active/`; sweeping a parallel feature's unresolved drafts into the wrong archive folder would be silent data displacement. Verified state first; parallel session had already self-archived. |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Compiled SuperPRD | `spec/compiled/SuperPRD_StereoPhase4_ZoneA_P.md` | created |
| MiniPRD-A (encode) | `spec/compiled/MiniPRD_ZoneA_P_Encode.md` | created |
| MiniPRD-B (generation) | `spec/compiled/MiniPRD_ZoneA_P_GenEnforce.md` | created |
| Archived Draft PRD | `spec/archive/20260703_115608_ZoneA_P_Encoding/Draft_PRD.md` | moved from `spec/active/` |
| Archived Red Team report | `spec/archive/20260703_115608_ZoneA_P_Encoding/RedTeam_Report_ZoneA_P_Encoding.md` | moved from `spec/active/` |
| Project memory index | (Claude Code project memory, outside repo) | updated stale SuperPRD pointer |

`spec/active/` is now empty (`.gitkeep` only).

## Results & Outcomes

All ten Red Team findings carry documented human decisions. The feature is fully specified and execution-ready: `SuperPRD_StereoPhase4_ZoneA_P.md` holds the pinned architecture (lone-pair convention, normative dummy-copy recipe, rdCIPLabeler end-to-end, re-embed-only enforcement) with a Resolved Trade-offs Log tracing every decision to its Red Team finding; the two MiniPRDs decompose it into ordered, executor-ready task lists (A: 14 tasks establishing the `_OIN_CIPCode_LP` contract; B: 8 tasks consuming it, blocked on A). The archival script confirmed success and flushed the active directory without touching the parallel Phase 3 effort.

## How to Reproduce

1. **Prerequisites:** repo at commit `7d90376` on `main`; `spec/active/` containing a Draft PRD and its Red Team report (produced by `/hyper-architect` then `/hyper-redteam` in separate context windows); `python3` available.
2. In a fresh context window (Red Team contamination rule), invoke `/hyper-resolve against spec/active/Draft_PRD.md and spec/active/RedTeam_Report_ZoneA_P_Encoding.md`.
3. The agent reads `.agents/skills/hyper-resolve/SKILL.md`, extracts the report's consolidated blocking table, and walks the state machine: answer the forced trade-offs (max 2 per turn) — blockers first, then highs, then the grouped NFR-defaults question, then the Candidate Artifact routing question.
4. After the last decision, the agent reads `.agents/schemas/` templates plus an existing compiled SuperPRD for format, writes `SuperPRD_<Feature>.md` + `MiniPRD_*.md` to `spec/compiled/`, then runs `python3 .agents/scripts/archive_specs.py <Feature_Name>` and logs the returned absolute path.
5. **Gotchas:** the archive script flushes *every* file in `spec/active/` — if a parallel feature's drafts live there, stash them first (or confirm the parallel session already archived them, as happened here). Expected end state: `spec/active/` empty, `spec/compiled/` holding the new SuperPRD + MiniPRDs, archive folder timestamped under `spec/archive/`.

## Patterns & Lessons

- **Cross-convention comparisons are a spec smell:** the PRD's own §5.1 declared the two CIP views legitimately divergent, yet §3.1.4 compared them. When a document defines a convention boundary, audit every comparison in the spec against it.
- **"Reuse machinery from a parallel draft" is a hidden blocker:** specifying against another unresolved draft's not-yet-built utility (Phase-3 reflection) couples two features' schedules invisibly. The resolve step should either decouple (as here — re-embed needs no reflection) or make the dependency an explicit blocker.
- **Improper transforms are never local:** any mirror-based "fix one stereocenter" plan inverts every stereocenter in the object. Prefer regeneration-from-constraint (re-embed from the tag) over geometric surgery.
- **Flush-style scripts + parallel features:** HACF's single `spec/active/` directory assumes one feature in flight; with two parallel drafts, the archival flush becomes destructive. Check directory contents immediately before flushing — state can change between reads when parallel sessions run (Phase 3 self-archived mid-session here).
- **Environment note:** `python` is not on PATH on this machine; use `python3` for `.agents/scripts/`.
