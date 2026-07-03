# Process Document: Stereo Roadmap Phase 1 — Winding Plumbing Spec Chain

**Generated:** 2026-07-03T08:13:48
**Session Focus:** Author and harden the HACF spec chain (Draft PRD → Red Team → compiled SuperPRD + MiniPRD) for ROADMAP-stereo.md Phase 1, "Preserve the Signal" — making eta-ligand winding direction survive the OIN→XYZ parse path. No production code was written; the deliverable is an execution-ready spec.

## Problem Statement

The OIN→XYZ generation pipeline silently discards eta-ligand **winding direction** (`{n>}` = clockwise heading, `{n<}` = counter-clockwise heading). Phase 0 diagnostics (TASK-10) proved this empirically: flipping ferrocene's ring-0 marker from `{0>}` to `{0<}` produces byte-identical 3D output. The signal is destroyed at the very first parse step — `SLOT_REGEX` in `src/oinsmiles/oin/inline.py:44` matches the `>`/`<` suffix but never captures it. Phase 1's job is plumbing only: capture the suffix and carry it as far as `ParsedOIN` so a later phase can act on it. No geometry behavior changes this phase.

## Starting State

- Git HEAD: `7edae02` ("Add Phase 0 stereo diagnostics; correct roadmap fixture assumptions"). Test suite green (`discover tests` and `discover tests/unit` both OK; one intentional `expectedFailure` = `test_haptic_face_winding`, the Phase 3 gate).
- The multi-session v3.7 + stereo effort is driven from `spec/worklog/NOTES.md`, with append-only decisions D-1…D-5 and a phase roadmap in `spec/worklog/ROADMAP-stereo.md`.
- Phase 0 (TASK-10) was DONE and had *revised the plan*: winding loss is the real live gap; ligand carbon `@/@@` already survives the template path; genuine P/N-center coverage needs a new fixture (deferred to Phase 2).
- `spec/active/` was empty — no in-progress PRD. The roadmap mandated that stereo phases "graduate" from roadmap direction to a full HACF MiniPRD (`/hyper-architect` → `/hyper-redteam` → `/hyper-resolve`) before any implementation.

## Approach & Methodology

Spec-driven, adversarial-loop HACF workflow, orchestrated under `/remote-control OIN-2`. The sequencing was the canonical three-phase chain: an opinionated interview to produce a Draft PRD (`/hyper-architect`), an adversarial review to find gaps (`/hyper-redteam`), and a resolution/compilation pass that dispositions every finding and emits execution-ready specs (`/hyper-resolve`). The rationale for the full chain (rather than implementing straight from the roadmap) is a standing roadmap rule: winding touches a load-bearing data contract consumed in multiple places, so the consumer map and the back-compat story needed adversarial verification before code.

## Steps Taken

1. **Read the context first** (`/hyper-architect`, codebase-first rule). Read `spec/worklog/NOTES.md`, `ROADMAP-stereo.md`, and — before asking the user anything — the actual producer and every consumer of the data structure in question. Grep for `parse_inline_string` / `vector_data` surfaced exactly four sites: the producer (`oin/inline.py:353`), a 3-way unpack (`generation/oin_parser.py:485`), a **2-way** unpack (`oin/parser.py:34`, already broken against today's 3-tuples — the dead TD-003 `SMILESToXYZ` path), and exact-equality asserts (`tests/unit/test_inline.py:21,30`).

2. **Surfaced the one genuine fork via AskUserQuestion.** The scope was precise enough that only the flagged main risk needed a human decision: how to change the `vector_data` element type. Options weighed — plain 4-tuple (positional, undocumented 4th slot), sidecar dict (second channel, desync risk), or a `typing.NamedTuple`. **User selected the `NamedTuple`** (`SlotAssignment(lig_rank, atom_idx, slot, winding=None)`): named `.winding` access, self-documenting, and still a real tuple so index reads `[0..2]` stay valid.

3. **Resolved the remaining choices as architect** (not user-blocking): store winding as the literal `'>'`/`'<'`/`None` char (no premature enum); scope threading to the inline path only (legacy V2.4 sidecar `w:` stays `None`); treat `_build_connected_smiles` as verify-unchanged; and hold `test_haptic_face_winding` at `expectedFailure` (it is Phase 3's acceptance gate, not Phase 1's).

4. **Wrote `spec/active/Draft_PRD.md`** using the SuperPRD structure, with a §5.1 consumer map as the centerpiece and a §9 note directing the Red Team at the three highest-risk areas (consumer-map completeness, `NamedTuple` pickling across the molassembler `ProcessPoolExecutor` boundary, and guarding against winding leaking into placement). Appended a Log entry to `NOTES.md`.

5. **Ran `/hyper-redteam`** on the Draft. Verdict: "unusually disciplined"; §5.1 map CONFIRMED complete against a repo-wide grep. Five findings, two of which changed acceptance criteria (see Key Decisions). Report saved alongside the draft.

6. **Ran `/hyper-resolve`** — dispositioned all five findings and compiled the output. Flushed `spec/active/` via `archive_specs.py`, moving the Draft PRD and Red Team report to `spec/archive/20260703_080740_StereoPhase1_Winding/`, and emitted the compiled specs to `spec/compiled/`. Updated `NOTES.md` session state to "COMPILED — ready for /hyper-execute."

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| `SlotAssignment` as a `typing.NamedTuple` with `winding=None` | Plain 4-tuple; sidecar winding dict alongside the 3-tuple | Named access + still a real tuple; a missed positional-index reader stays correct, only positional-*unpack* fails fast (intended) |
| Winding preserved on a new `ParsedOIN.winding_by_slot` dict, populated on **all** paths | Copy winding onto `OINVector` at `oin_parser.py:485` (the Draft's plan) | **RT #1 (material gap):** that copy sits inside `if tmpl_vectors is not None:`, so template-less geometries (`NON`, template-less eta — the exact haptic/eta family Phase 3 targets) got `vectors=[]` and lost winding. Ferrocene passed only because `LIN` is a template key. A `ParsedOIN`-level channel keeps the adapter's `vectors` iteration byte-inert *structurally*, not just by test. |
| Parse regex `\{(\d+)([><^])?\}`, normalize `^`→`>` on capture | Draft's `\{(\d+)([><])?\}` (only `>`/`<`) | **RT #2:** the suffix is the *heading* marker; the generate side already normalizes `^`→`>` (`oin/inline.py:245`). Missing `^` would silently drop a legal marker. |
| Store literal `>`/`<`/`None`; no enum/bool | Normalize to a direction enum or boolean | Faithful to source and to the Phase-0 diagnostic strings; avoids premature abstraction |
| Add a pre/post XYZ byte-diff harness as a hard §8 gate | Infer inertness from a green suite (Draft's approach) | **RT #5:** "no geometry change" must be *proven*, not assumed from unrelated green tests |
| Drop the `NamedTuple` pickling stress item | Keep it as a Red Team focus area | **RT #4:** the `ProcessPoolExecutor` worker boundary (`molassembler_adapter.py:2220`) carries a primitives-only dict; `OINVector`/`SlotAssignment` never pickle. Non-issue. |
| Hold `test_haptic_face_winding` at `expectedFailure` | Wire winding into placement to flip it green | Phase 1 is plumbing; that test is Phase 3's acceptance gate. Flipping it now would be an out-of-scope behavior leak. |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Draft PRD | `spec/active/Draft_PRD.md` → archived to `spec/archive/20260703_080740_StereoPhase1_Winding/Draft_PRD.md` | created, then archived |
| Red Team report | `spec/archive/20260703_080740_StereoPhase1_Winding/RedTeam_Report.md` | created |
| Compiled SuperPRD | `spec/compiled/SuperPRD_StereoPhase1_Winding.md` (v1.0.0, confidence 10/10) | created |
| Compiled MiniPRD | `spec/compiled/MiniPRD_WindingPlumbing_Phase1.md` (12 tasks, 5 verification tests) | created |
| Worklog | `spec/worklog/NOTES.md` | updated (architect Log entry + redteam/resolve Log entry + session-state advanced to "ready for /hyper-execute") |

No files under `src/` were modified in this session.

## Results & Outcomes

An execution-ready spec chain for Phase 1 now exists in `spec/compiled/`: `SuperPRD_StereoPhase1_Winding.md` and `MiniPRD_WindingPlumbing_Phase1.md` (12 implementation tasks, 5 verification tests, confidence 10/10). The Red Team turned up one genuinely material gap the Draft had missed — winding would have been *conditionally* dropped at `ParsedOIN` for exactly the template-less eta geometries Phase 3 needs — which is now fixed at the design level by routing winding through a `ParsedOIN.winding_by_slot` dict populated on all paths, outside the template gate. The suite remains green and `test_haptic_face_winding` remains a deliberate `expectedFailure` (the Phase 3 gate). The next action is `/hyper-execute MiniPRD_WindingPlumbing_Phase1.md`.

## How to Reproduce

Prerequisite state: git at `7edae02` (or later with Phase 0 diagnostics present), `uv` toolchain installed, suite green (`uv run python -m unittest discover tests` and `.../discover tests/unit` both OK).

1. **`/hyper-architect`** with the Phase 1 scope (SLOT_REGEX winding capture; `vector_data` element gains a winding field default `None`; thread to `ParsedOIN`; no placement change). Read `NOTES.md` + `ROADMAP-stereo.md` first, then grep every consumer of `parse_inline_string`/`vector_data` before asking anything. Expect one interview question (the tuple-contract fork). → produces `spec/active/Draft_PRD.md`.
2. **Start a new conversation** (per CLAUDE.md, the Red Team must not inherit the Architect's history) and run **`/hyper-redteam`** on `spec/active/Draft_PRD.md`. → produces a Red Team report; expect it to confirm the consumer map and probe the `ParsedOIN` preservation path.
3. **`/hyper-resolve`** to disposition every finding and compile. → emits `SuperPRD_*` + `MiniPRD_*` to `spec/compiled/` and archives the draft + report to `spec/archive/<timestamp>_StereoPhase1_Winding/`.
4. **`/hyper-execute MiniPRD_WindingPlumbing_Phase1.md`** to implement (not done in this session).

Gotchas / order-dependencies: run architect and redteam in **separate conversations** (adversarial isolation). Under `/remote-control`, the active PRD is archived by `/hyper-resolve`, so read compiled specs from `spec/compiled/` and archived source-of-record from `spec/archive/`, not `spec/active/` (which is flushed empty). Verify with `find` if a `spec/active/` artifact appears "missing" — it was likely archived, not lost.

## Patterns & Lessons

- **Codebase-first interviewing collapses the question count.** Enumerating the four `vector_data` consumers before the interview meant only the one true fork (tuple representation) needed the user; everything else was a defensible architect default.
- **Test-passing is not proof of design correctness.** The Draft's per-`OINVector` winding copy passed the ferrocene case purely because `LIN` happens to be a template key. The Red Team's value was catching the *conditional* silent drop that no current fixture would expose — the fix (a `ParsedOIN`-level channel outside the template gate) is structurally inert rather than merely tested-inert.
- **The heading marker is a three-symbol alphabet (`>`/`<`/`^`), not two.** A regex derived from one direction of the round-trip missed the normalization the generate side already does. Always reconcile a new parse regex against the corresponding serializer.
- **Guard against scope leak explicitly.** Making "hold `expectedFailure`" and "prove XYZ byte-identity" *acceptance criteria* prevents a well-meaning implementer from wiring winding into placement to turn the Phase 3 gate green early.
- **`/remote-control` orchestration archives as it compiles.** A file written to `spec/active/` in an earlier turn can legitimately be gone by a later turn — check `spec/archive/` before treating it as a lost write.
