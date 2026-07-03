# Process Document: Stereo Phase 3 — Haptic Face Control (Architect session)

**Generated:** 2026-07-03T11:23:29
**Session Focus:** Running `/hyper-architect` to produce the Draft PRD for ROADMAP-stereo.md Phase 3 (haptic face control on the OIN→XYZ generation path), and resolving the roadmap's open design question about V3.6 winding semantics.

## Problem Statement

The OIN→XYZ generation direction silently ignores eta-ring winding: flipping a ring's `{n>}`/`{n<}` marker in an input OIN produces byte-identical 3D output, so the prochiral face of a coordinated Cp/arene is uncontrollable. Phase 1 already plumbed the winding marker as far as `ParsedOIN`; this session's goal was to architect Phase 3 — the placement-side correction that finally *consumes* that winding — and to resolve the roadmap's standing open question: does V3.6 winding pin the prochiral face for substituted rings, or only ring direction?

## Starting State

- Git HEAD: `7d90376` (Phase 4 design brief + parallel Phase 3/4 plan). Branch `main`, ahead of origin.
- Stereo Phase 0 (diagnostics) and Phase 1 (winding plumbing) DONE and committed. Winding reaches `ParsedOIN.winding_by_slot` and per-vector `OINVector.winding`, but nothing on the generation path reads it.
- `tests/unit/test_stereo_roundtrip_diagnostics.py::test_haptic_face_winding` is an `@unittest.expectedFailure` documenting the gap (ferrocene, `OUT A == OUT B`).
- The eta placement code (`_stitch_eta_fragment`, `_stitch_multi_eta_fragment` in `src/oinsmiles/generation/molassembler_adapter.py`) does centroid-plane SVD alignment that fixes the ring normal but leaves circulation sense arbitrary.
- A parallel Phase 4 (Zone-A P stereocenter) session was already running and had written `spec/active/Draft_PRD.md`.
- The working tree carried a large unrelated ruff-format churn owned by another session (per NOTES.md ⚠️) — to be left untouched.

## Approach & Methodology

Spec-driven, single HACF phase: `/hyper-architect` only. The Architect role mandates codebase-first exploration — resolve everything derivable from source before interviewing, then ask one question at a time with a recommended default. The sequencing was depth-first: core mutation → data/blast-radius → personas → novel frontier → draft. The distinguishing feature of this session was that codebase analysis surfaced a *geometric impossibility* that contradicted the task brief, so the interview became a series of "here is what the code/math forces; do you accept the correction?" confirmations rather than open elicitation.

## Steps Taken

1. Read `NOTES.md`, `ROADMAP-stereo.md`, and the diagnostic test to load the multi-session state and the exact acceptance target. Observed the roadmap flags an open question about winding semantics and suggests a methylcyclopentadienyl fixture may be needed.
2. Invoked `/hyper-architect`; read its `SKILL.md` (relentless one-question-at-a-time interviewer, codebase-first).
3. Explored source before asking anything: `molassembler_adapter.py` (`_stitch_eta_fragment` `:584`, `_stitch_multi_eta_fragment` `:151`, caller `_template_generate` `:1034`, post-placement optimiser `:1170`); `oin_aligner.py::_determine_winding` `:634`; `oin_parser.py` (`OINVector.winding`, `ParsedOIN.winding_by_slot`). Established that winding is available two ways and that the stitch functions currently receive no winding.
4. **Phase 1 (core mutation):** Reasoned through the geometry of the roadmap's literal instruction ("mirror the fragment across the ring plane") and found it is a *no-op for planar rings* — reflection across a ring's own plane leaves in-plane atoms fixed, so it cannot reverse the measured circulation. Recommended re-framing the core operation as a face-swap that reverses circulation about the metal→centroid axis. User accepted.
5. **Phase 2, Q1 (winding channel):** Recommended consuming per-vector `OINVector.winding` rather than `winding_by_slot` (the brief's named source), because the integer slot is not in hand where eta groups are keyed by rounded vector-direction tuple. User accepted.
6. **Phase 2, Q2 (multi-eta scope):** Determined a bridged ansa-metallocene cannot have one ring corrected without tearing the Si bridge; recommended full correction for single-eta, coherent whole-fragment-only for bridged. User accepted and added that independent bridged-ring flips should later use tether-dihedral rotation.
7. **Phase 2, Q3 (the pivotal finding):** Proved that for **unsubstituted ferrocene the winding is not a geometric observable** — a regular Cp is a set-invariant under the face-swap, so the swapped structure is geometrically identical and re-encodes identically (matching the existing `OUT A == OUT B` diagnostic). This makes `test_haptic_face_winding` unsatisfiable by any geometry-only correction and **resolves the roadmap's open question**: winding pins a physical face only for *substituted* rings. Recommended re-targeting acceptance to a desymmetrized fixture and demoting the ferrocene test to a documented skip. User agreed and supplied a superior fixture, `tests/integration/Ferrocene-halide-face.xyz` (each Cp bears five distinct substituents H/OH/Cl/Br/I).
8. **Phase 2, Q4 (operator):** Recommended a proper 180° in-plane rotation (det +1) over a reflection (det −1), because the reflection inverts pendant-substituent chirality and would collide with Phase 4. User accepted, noting it is the free-ligand degenerate case of a tether-dihedral flip.
9. **Phase 2, Q5 (convention parity):** Recommended extracting a shared `signed_circulation` helper as the single source of truth for the winding sign, called by both `_determine_winding` and the new correction, to eliminate the sign-inversion bug class. User accepted.
10. **Phase 3 (personas):** Derived N/A (pure in-process library transform). Advanced without a question.
11. **Phase 4 (novel frontier):** Ran `XYZToSMILES().convert()` on the fixture to ground the discussion empirically — confirmed it encodes with live winding markers on both rings. Recommended the Candidate-Artifact golden plus a **faithful per-ring round-trip** acceptance assertion (per-ring regex, robust to heading-atom identity) with a documented weaker geometric fallback. User accepted.
12. **Phase 5 (draft):** Read the SuperPRD schema. Attempted to write `spec/active/Draft_PRD.md` and discovered it was **already occupied by the parallel Phase 4 session**. Did not clobber; wrote the Phase-3 draft to `spec/active/Draft_PRD_StereoPhase3_HapticFace.md` instead and flagged the collision. Appended a Log entry to `NOTES.md`.

## Key Decisions & Rationale

| Decision | Alternatives Considered | Reason Chosen |
|---|---|---|
| Consume `OINVector.winding` at placement | `ParsedOIN.winding_by_slot` (per brief) | Integer slot unavailable where eta groups are keyed by direction tuple; heading atom's `.winding` already carries the marker |
| Proper 180° in-plane rotation correction | Reflect-across-ring-plane (roadmap wording); reflection-through-axis-plane | Only operator that reverses circulation, swaps the face, AND preserves pendant chirality (det +1) — forward-compatible with Phase 4 |
| Extract shared `signed_circulation` helper | Duplicate the sign math in the adapter | Single source of truth; kills encode/generate sign-drift bug class |
| Re-target acceptance to substituted fixture; skip ferrocene test | Force the ferrocene test to pass | Ferrocene winding is not a geometric observable — proven unsatisfiable; faking it would encode a non-physical distinction |
| Single whole-fragment correction for bridged multi-eta | Independent per-ring correction | Can't reflect one ring without tearing the Si bridge; independent control deferred to tether-dihedral rotation |
| Faithful per-ring round-trip assertion | Weaker "outputs differ" gate | Convention parity makes the strong assertion achievable; per-ring regex handles heading-atom variance |
| Save draft under distinct filename | Overwrite `spec/active/Draft_PRD.md` | Slot owned by the concurrent Phase 4 session — must not clobber |

## Artifacts Created / Modified

| Artifact | Path | Change |
|---|---|---|
| Phase 3 Draft PRD | `spec/active/Draft_PRD_StereoPhase3_HapticFace.md` | created |
| Worklog Log entry | `spec/worklog/NOTES.md` | updated (appended Phase 3 architect entry) |
| This process document | `spec/process/process_20260703_112329_stereo-phase3-haptic-face-architect.md` | created |
| Acceptance fixture (user-provided) | `tests/integration/Ferrocene-halide-face.xyz` | pre-existing this session (referenced) |

No `src/` code was modified. No git commit was made.

## Results & Outcomes

- A complete, RedTeam-ready Draft PRD for Phase 3 exists, with all five interview branches resolved and a 9/10 confidence score.
- The roadmap's long-standing open question is **resolved**: V3.6 winding pins a physically observable prochiral face only for substituted rings.
- The acceptance strategy is corrected from an unsatisfiable target (ferrocene hard pass) to an achievable one (faithful per-ring round-trip on the desymmetrized halide fixture, verified to encode as `[Fe_LIN].Oc{0<}1[cH]{0}c{0}(Cl)c{0}(Br)c{0}1I.Oc{1}1[cH]{1}c{1}(I)c{1<}(Br)c{1}1Cl`).
- Two candidate MiniPRDs are scoped (shared helper; face correction) for `/hyper-resolve` to compile.

## How to Reproduce

1. **Prerequisite state:** repo at/after `7d90376` on `main`; Stereo Phase 1 committed (`6820d3a`); `uv sync` done; `tests/integration/Ferrocene-halide-face.xyz` present.
2. Read `spec/worklog/NOTES.md` and `spec/worklog/ROADMAP-stereo.md` to load multi-session state.
3. Invoke `/hyper-architect` with the Phase 3 scope. Expect the Architect to explore `molassembler_adapter.py`, `oin_aligner.py`, and `oin_parser.py` before asking anything.
4. Answer the five design questions (winding channel, multi-eta scope, fixture re-target, operator, convention parity). The winding-observability proof should independently reproduce: a regular Cp is set-invariant under the face-swap.
5. Empirically verify the fixture: `uv run python -c "from oinsmiles import XYZToSMILES; print(XYZToSMILES().convert('tests/integration/Ferrocene-halide-face.xyz'))"` → OIN with winding markers on both rings.
6. Expect output `spec/active/Draft_PRD_StereoPhase3_HapticFace.md`. **Gotcha:** if a parallel session occupies `spec/active/Draft_PRD.md`, save under a distinct name — do not overwrite.
7. **Next phase:** new conversation → `/hyper-redteam` against the Phase-3 draft file (not the default `Draft_PRD.md`).

## Patterns & Lessons

- **Codebase-first can overturn the brief.** Geometric/mathematical analysis of the placement code revealed the assigned acceptance test was physically unsatisfiable. The Architect's job here was less requirement elicitation than surfacing a hard truth and getting the deliverable re-scoped — with explicit user sign-off on the deviation.
- **"Mirror across the ring plane" was imprecise.** The correct operator is a proper 180° in-plane rotation, not a reflection; the distinction matters for chirality preservation and interacts with a downstream phase. Always work the actual transform math, not the prose.
- **Symmetry kills observability.** A stereo signal that a symmetric object cannot express should be tested on a desymmetrized fixture — otherwise the test measures nothing.
- **Parallel-session hygiene.** With concurrent Phase 3/4 sessions sharing `spec/active/`, always check the target before writing a canonically-named artifact and flag collisions rather than clobbering.
- **Choose a shared helper when two modules must agree on a sign/convention** — it converts a whole class of silent-inversion bugs into a compile/grep-visible single definition.
