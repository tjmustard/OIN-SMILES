# Process Document: Stereo Roadmap Phase 1 — Red Team Resolution & Spec Compilation

**Generated:** 2026-07-03T08:14:55-07:00
**Session Focus:** Mediating the Red Team report for Stereo Phase 1 ("Preserve the Signal" — winding plumbing) and compiling the final SuperPRD + MiniPRD via `/hyper-resolve`.

## Problem Statement

The OIN→XYZ generation pipeline silently discards eta-ligand **winding direction** (`{n>}` = CW heading, `{n<}` = CCW heading), because the parse regex matches but never captures the direction suffix. A Draft PRD and adversarial Red Team review already existed for the fix; this session's job was to mediate the Red Team's findings into definitive architectural decisions and compile the executable specification.

## Starting State

At commit `7edae02` (Phase 0 stereo diagnostics shipped). The relevant inputs were in `spec/active/`:
- `Draft_PRD.md` (v0.1.0) — a tightly scoped "plumbing only" design produced by `/hyper-architect`, resolving the core fork to a `SlotAssignment` `NamedTuple`.
- `RedTeam_Report.md` — an adversarial review produced by `/hyper-redteam`, verdict "unusually disciplined," with the §5.1 consumer map CONFIRMED complete against a repo-wide grep. It raised five findings, two of which changed acceptance criteria.

No compiled Phase 1 spec existed yet. The generation-side winding gap was empirically proven by `test_haptic_face_winding` (Phase 0), which stays `@unittest.expectedFailure` until Phase 3.

## Approach & Methodology

Spec-driven, following the HACF Phase 1 chain: this session was the `/hyper-resolve` node (`/hyper-architect` → `/hyper-redteam` → **`/hyper-resolve`** → ready for `/hyper-execute`). Methodology: read both active artifacts, triage the Red Team findings by severity, drive the two acceptance-changing decisions through forced-choice `AskUserQuestion` prompts (max 2 per turn per skill rules), apply sensible defaults to the three hardening items, then compile and archive.

## Steps Taken

1. **Read `RedTeam_Report.md` + `Draft_PRD.md`** to establish the finding set. Observed the Red Team had already cleared two of three original stress items (pickle boundary, §5.1 completeness) and surfaced one net-new material gap not anticipated by the Draft's Execution Note.

2. **Triaged into two severity tiers.** The acceptance-changing findings were #1 (conditional winding drop) and #2 (`^` alphabet); #3/#4/#5 were wording/hardening. Presented #1 and #2 as forced trade-offs via `AskUserQuestion`.
   - **#1 outcome:** user chose "build vectors on NON path" (deliver the guarantee for template-less geometries).
   - **#2 outcome:** user chose three-symbol alphabet — accept `^`, normalize `^`→`>`.

3. **Detected a collision between the #1 fix and the Draft's own inertness constraint.** Building `OINVector`s on the `tmpl_vectors is None` path would make the adapter's `_build_connected_smiles` see vectors where it previously saw `[]` — breaking §6 "zero diff in generated XYZ." Raised this as a follow-up forced choice rather than silently compiling a contradiction.
   - **Outcome:** user chose "separate winding channel" — a `ParsedOIN.winding_by_slot` dict populated on all paths, leaving `NON`-path `vectors` emission untouched. This preserves the §1.2 guarantee (winding at `ParsedOIN` for all geometries) **and** structural inertness. This superseded the literal "placeholder OINVector" mechanism from step 2 while keeping its promise.

4. **Batched the three hardening findings** (#3 wording, #4 delete pickle stress-item, #5 byte-diff harness) into one "approve defaults" `AskUserQuestion`. User approved all.

5. **Compiled the SuperPRD** to `spec/compiled/SuperPRD_StereoPhase1_Winding.md` (v1.0.0), including a §9 Red Team disposition log mapping each finding to its decision.

6. **Compiled the MiniPRD** to `spec/compiled/MiniPRD_WindingPlumbing_Phase1.md` — 12 atomic tasks, 5 verification tests, confidence 10/10.

7. **Ran the archival script** `python3 .agents/scripts/archive_specs.py StereoPhase1_Winding`, flushing `spec/active/` → `spec/archive/20260703_080740_StereoPhase1_Winding/` (Draft + Red Team preserved).

8. **Updated `spec/worklog/NOTES.md`** — Phase 1 row set to "COMPILED — ready for /hyper-execute", "Next up" pointer updated, and a full Log entry recording all five dispositions.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| Universal `ParsedOIN.winding_by_slot` dict, populated on all paths outside the template gate | (a) Emit placeholder `OINVector`s on the NON path; (b) document winding as template-geometry-only | (a) perturbs the adapter's `vectors` iteration → breaks inertness; (b) abandons the exact eta/haptic family Phase 3 targets. The dict delivers the guarantee AND keeps `vectors` emission byte-inert. Also immune to the `slot_idx` overflow drop (keyed by slot). |
| Three-symbol winding alphabet: regex `([><^])?`, normalize `^`→`>` on capture | Two-symbol `([><])?` (reject `^`) | The generate side already normalizes `^`→`>` (`oin/inline.py:245`); a legacy `{0^}` would otherwise leak a literal `^` into output SMILES and yield `winding is None`. Symmetric handling closes the gap. |
| Byte-diff XYZ harness as a hard gate on inertness | Infer inertness transitively from a green suite | A green `unittest discover` does not prove byte-identical XYZ if no test asserts on exact coordinates. The harness converts §8's strongest claim from assertion to evidence. |
| Delete the pickle / `ProcessPoolExecutor` stress item | Keep it as a residual risk | Red Team verified the worker boundary carries a primitives-only dict; `OINVector`/`SlotAssignment` never cross it. Genuine non-issue. |
| RISK-1/US-003 reworded to "closed enumerated set + index-safe/unpack-fail-fast" | Keep "graceful positional degradation" | The NamedTuple protects `sa[0..2]` index reads but a 4-field arity makes `a,b,c = sa` fail fast (desirable). "Graceful degradation" mis-framed the real safety, which is the grep-closed consumer set. |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| SuperPRD (Phase 1) | `spec/compiled/SuperPRD_StereoPhase1_Winding.md` | created |
| MiniPRD (Phase 1) | `spec/compiled/MiniPRD_WindingPlumbing_Phase1.md` | created |
| Draft PRD | `spec/active/Draft_PRD.md` → `spec/archive/20260703_080740_StereoPhase1_Winding/` | archived |
| Red Team Report | `spec/active/RedTeam_Report.md` → `spec/archive/20260703_080740_StereoPhase1_Winding/` | archived |
| Worklog | `spec/worklog/NOTES.md` | updated (status table, next-up, Log entry) |

## Results & Outcomes

Phase 1 now has a compiled, execution-ready specification with all five Red Team findings dispositioned (confidence 10/10). The material correctness gap — winding silently dropped for template-less (`NON`/eta) geometries, the family Phase 3 depends on — is resolved by design (`winding_by_slot`) rather than deferred. `spec/active/` is flushed to prevent context collapse. The next actionable step is `/hyper-execute` against `MiniPRD_WindingPlumbing_Phase1.md`. No `src/` code was modified this session.

## How to Reproduce

Prerequisite: repo at a commit where `spec/active/` holds a `Draft_PRD.md` and a `RedTeam_Report.md` for the feature (i.e., after `/hyper-architect` and `/hyper-redteam`).

1. Invoke `/hyper-resolve`. It reads `.agents/skills/hyper-resolve/SKILL.md`, then `spec/active/RedTeam_Report.md` and `Draft_PRD.md`.
2. Triage: the agent presents the highest-severity, acceptance-changing findings first as forced-choice `AskUserQuestion` prompts (≤2 per turn). Answer them.
3. Resolve any decision that collides with an existing constraint (here: NON-path vectors vs. inertness) before compiling — do not silently compile a contradiction.
4. Approve or modify standard defaults for the remaining hardening findings.
5. The agent writes `spec/compiled/SuperPRD_<Feature>.md` and `spec/compiled/MiniPRD_<Module>.md`.
6. The agent runs `python3 .agents/scripts/archive_specs.py <Feature_Name>` — note `python3` (the environment has no bare `python`). Output prints the absolute archive path; `spec/active/` is left with only `.gitkeep`.
7. Update `spec/worklog/NOTES.md` with a Log entry.

Gotcha: the archival script is order-dependent — run it **after** the compiled specs are written, since it flushes everything in `spec/active/` regardless of content.

## Patterns & Lessons

- **Compose, don't overwrite, user choices.** The step-2 "build vectors on NON path" answer and the step-3 "separate channel" answer were reconciled into a single coherent design (guarantee via dict + inert vectors) rather than treating the later answer as a contradiction of the earlier. Surface the collision explicitly and let the user pick the reconciliation.
- **A green test suite is not a byte-diff.** When a constraint is "zero output change," add an explicit pre/post byte-diff harness; do not rely on an unrelated green suite as evidence.
- **The most dangerous finding is the one the author's own Execution Note didn't anticipate** — here, the conditional winding drop hidden behind a lucky templated fixture (ferrocene/`LIN`). Prioritize findings that change acceptance criteria over wording fixes.
- **Environment note:** this repo's scripts run under `python3`, not `python` (no bare `python` on PATH).
