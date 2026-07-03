# OIN-SMILES Worklog — Session Notes

Persistent state for the multi-session v3.7 + test-recovery + stereo effort.
Any session (human, Claude, or subagent) should read this file first and append
a Log entry before ending.

**Process:** Each `TASK-NN-*.md` file in this directory is a self-contained spec
implementable by a cheap model (Haiku/Sonnet/Opus) with zero other context.
Implementer instructions: read the task file, follow it exactly, run the
acceptance commands, then update the task's `Status:` line and add a Log entry
here. Stereo feature work (Phases 1–4 in `ROADMAP-stereo.md`) graduates to full
HACF MiniPRDs (`/hyper-architect` → `/hyper-redteam` → `/hyper-resolve`) when
its phase begins — do not implement those from the roadmap alone.

---

## Session state

Test baseline (2026-07-02, before any fixes):
- `uv run python -m unittest discover tests` → 60 collected: 7 failures + 5 errors
- `uv run python -m unittest discover tests/unit` → 52 run: 50 pass, 1 fail, 1 skip
- Integration dir (`tests/integration`) aborts on discovery (HACF framework test
  requires Python 3.11) — **out of scope by decision**

**Current status (2026-07-03, after TASK-01/02/03/04/10): GREEN + v3.7 shipped.**
- `uv run python -m unittest discover tests` → 55 run, OK
- `uv run python -m unittest discover tests/unit` → 55 run, OK (skipped=1,
  expected failures=1) — includes the 3 new stereo diagnostics
- Note: root discovery does not recurse into `tests/unit` (no `__init__.py`
  there — see AGENTS.md); run both commands to cover everything.
- Committed: 6950bff (tests green), d691a5f (v3.7). TASK-10 + roadmap update
  committed after.
- **Stereo Phase 1 (winding plumbing) is DONE** (2026-07-03): full HACF chain
  (`/hyper-architect` → `/hyper-redteam` → `/hyper-resolve` → `/hyper-execute`
  → `/hyper-audit`), plus a post-audit review fix (winding_by_slot clobbering
  guard — see Log). MiniPRD archived as
  `spec/archive/MiniPRD_WindingPlumbing_Phase1_AUDITED.md`. Suite: 63 unit
  tests OK (1 skip, 1 expected failure = the haptic diagnostic, which stays
  red until Phase 3 *uses* the winding).
- **Next up: Phase 3 (haptic face control) or the Phase-2 prerequisite fixture**
  (a genuine P/N-stereocenter complex — see roadmap Phase 0 + Fixtures).
  Phase 3 can start now that winding reaches `ParsedOIN.winding_by_slot`;
  needs its own MiniPRD via the HACF chain.

| Task | Status | Model tier | Depends on |
|---|---|---|---|
| TASK-01 delete stale root tests | DONE | Haiku | — |
| TASK-02 vertex-index expectations | DONE | Haiku | — |
| TASK-03 GeneratedStructure assert | DONE | Haiku | — |
| TASK-04 v3.7 descriptor-free metal token | DONE | Sonnet | TASK-01 |
| TASK-10 stereo diagnostic round-trips | DONE | Sonnet | TASK-04 |
| Stereo Phase 1 (winding plumbing) | DONE (executed + audited + review fix) | Sonnet/Opus | TASK-10 |
| Stereo Phases 2–4 | ROADMAP (needs MiniPRDs) | Sonnet/Opus | Phase 1 |

## Decisions (append-only)

- **D-1 (2026-07-02) — "v3.7" = descriptor-free metal token; the `@SP1`/`@OH10`
  segment is a bug, not a feature.** Git archaeology: true v3.6 (`68854a4`)
  emitted `[Pt_SPL]`; the `@desc` appeared at `0711d06` (molassembler
  integration) and the v4.0 design docs (`299a073` docs/PRD.md) specify
  descriptor-free. Root cause: stale `is_metal` local in
  `src/oinsmiles/utils/xyz2mol.py` — set in the first fragment loop (:868),
  never reassigned in the second loop (:910), so the metal fragment takes the
  ligand SMILES branch (:982-986) and RDKit's non-tetrahedral stereo extension
  (`CHI_SQUAREPLANAR` etc., set by `AssignAtomChiralTagsFromStructure` in
  `core/chirality.py:98`) serializes as `[Pt@SP1]`. The intended `[Pt]` branch
  (:988) is dead code. Removal loses nothing: cis/trans and fac/mer are fully
  encoded by slot ordering, and `parse_inline_string` discards `@desc` anyway
  (`oin/inline.py:338`).
- **D-2 (2026-07-02) — `vertex_indices` are 1-based fragment ranks; fix the
  tests, not `_extract_oin_constraints`.** Fragment rank 0 = the metal by
  invariant (`generation/engine.py:158`; ligands enumerate from `start=1`).
  Re-basing to 0 would make ligand entries resolve to the metal atom and trip
  the bond-to-self guard (`generation/oin_parser.py:180`).
- **D-3 (2026-07-02) — parsers stay tolerant of legacy `@desc` forever.**
  `METAL_REGEX` (`oin/inline.py:41`) already treats it as optional; keep it so
  old strings remain parseable. `tests/test_direct_parser_regex.py` feeds
  literal `@SP1` inputs and doubles as the back-compat coverage — keep those
  inputs as-is.
- **D-4 (2026-07-02) — generation currently DROPS stereo silently; fix
  information flow before evaluating alternative builders.** OIN→XYZ discards
  winding `{n>}`/`{n<}` (SLOT_REGEX doesn't capture it), never reads P/N CIP
  codes (`core/chirality.py` is XYZ→OIN only), and has no haptic-face control
  (only a ring-normal flip toward the metal). Molassembler is only the
  fallback path; primary is RDKit ETKDG + Kabsch templates. Decision on an
  alternative builder is DEFERRED until Phases 1–3 of `ROADMAP-stereo.md`
  produce data — the current builders have never been given the stereo signals.
- **D-5 (2026-07-02) — HACF framework tests excluded.** The failing
  `tests/integration` items (Python 3.11 gate, pytest-style fixtures) test the
  HACF toolchain, not this library.

## Log

### 2026-07-02 — Planning session (Fable)
- Explored: format history, full test-suite status, generation-side stereo
  handling, git archaeology of the `@desc` drift. Findings captured in D-1…D-5.
- Empirically verified (Plan agent, `sys.settrace` + monkeypatched pipeline
  runs): the `@desc` first appears as `sanitized_smiles='[Pt@SP1]'` at
  `xyz2mol.py:985`; simulated removal leaves slot ordering byte-identical and
  cis/trans, fac/mer, ligand `@/@@` all still distinguished.
- Created this worklog: TASK-01…04, TASK-10, ROADMAP-stereo.md.
- No source code was modified this session.

### 2026-07-03 — TASK-02 fix: vertex-index expectations (Haiku)
- Fixed 7 test expectations in `tests/test_direct_parser_regex.py` to match 1-based fragment rank convention
  - Line 20 (`test_platinum_square_planar_basic`): `[0, 1, 2, 3]` → `[1, 2, 3, 4]`
  - Line 58 (`test_no_shape_code`): `[0, 1]` → `[1, 2]`
  - Line 71 (`test_no_chiral_tag`): `[0, 1]` → `[1, 2]`
  - Line 42 (`test_iron_linear_with_heading_markers`): `[0, 1]` → `[1, 2]`
  - Line 112 (`test_heading_marker_clockwise`): `[0, 1]` → `[1, 2]`
  - Line 123 (`test_heading_marker_counterclockwise`): `[0, 1]` → `[1, 2]`
  - Line 132 (`test_multiple_vertex_indices`): `[0, 1, 2, 3, 4, 5]` → `[1, 2, 3, 4, 5, 6]`
  - Added comment documenting convention: `vertex_indices values are fragment ranks (metal = rank 0, ligands start at 1); list position = slot`
- Acceptance results: `uv run python -m unittest tests.test_direct_parser_regex -v` → **Ran 14 tests ... OK** ✓
- Task complete; no source modifications.

### 2026-07-03 — TASK-03 fix: GeneratedStructure assertion (Haiku)
- Fixed stale test assertion in `tests/unit/test_molassembler_adapter.py:133`
  - Changed `self.assertEqual(result, expected_xyz)` to `self.assertEqual(result.xyz, expected_xyz)`
  - Updated docstring to clarify return type is `GeneratedStructure` whose `.xyz` field holds the XYZ block
- Acceptance results:
  - `uv run python -m unittest tests.unit.test_molassembler_adapter -v`: OK (6 tests)
  - `uv run python -m unittest discover tests/unit`: 52 run, 0 failures, 1 skip ✓
- Task complete; no source modifications.

### 2026-07-03 — TASK-01 complete: Delete stale root-level tests (Haiku)
- **Orphan check verified:** `grep -rn "test_helpers" tests/ --include="*.py" | grep -v tests/unit` confirmed `tests/test_helpers.py` imported only by the 5 files being deleted.
- **Files deleted via `git rm`:**
  - `tests/test_axial_chiral.py`
  - `tests/test_binap_stability.py`
  - `tests/test_chiral_n.py`
  - `tests/test_chiral_p.py`
  - `tests/test_regression_stability.py`
  - `tests/test_helpers.py`
- **Acceptance results:**
  - `uv run python -m unittest discover tests | tail -5`: Ran 55 tests, FAILED (failures=7) — **0 import errors** ✓
  - `uv run python -m unittest discover tests/unit | tail -3`: Ran 52 tests, FAILED (failures=1, skipped=1) — **baseline unchanged** ✓
- No contradictions found. Task status set to DONE.

### 2026-07-03 — Reconciliation after parallel TASK-01/02/03 (Fable)
- TASK-01's acceptance run overlapped with the parallel TASK-02/03 agents, so
  the failure counts in its entry above are stale snapshots, not regressions.
- Combined verified state: `discover tests` → **55 run, OK**;
  `discover tests/unit` → **52 run, OK (skipped=1)**. Suite fully green.
- Remaining known-red: none (axial-chiral skip is intentional; integration dir
  excluded per D-5).

### 2026-07-03 — TASK-04 complete: v3.7 descriptor-free metal token (Sonnet)
- **Source fix (1 line):** `src/oinsmiles/utils/xyz2mol.py` — in the second
  per-fragment loop (`for i, item in enumerate(fragments_data):`), added
  `is_metal = item['is_metal']` immediately after `indices = item['indices']`.
  This stops the metal fragment from taking the stale `is_metal=False` ligand
  branch, so it now flows through the intended `sanitized_smiles = f"[{...}]"`
  path → `[Pt_SPL]` instead of `[Pt@SP1_SPL]`.
- **Confirmation step:** ran `discover tests/unit` immediately after the source
  fix (before touching expectations) — got exactly 7 failures (`test_n_stability`,
  `test_p_stability`, `test_cis_ptcl2en`, `test_cisplatin`, `test_fac_irppy3`,
  `test_mer_irppy3`, `test_transplatin`), every diff being *only* the removed
  `@desc` segment (verified via diff inspection) — no slot-order or other
  changes, so no STOP condition triggered.
- **Test expectations updated** (transformation rule: strip `\[([A-Z][a-z]?)@[A-Z]+[0-9]+_` → `[\1_`):
  - `tests/unit/test_regression_stability.py` — 5 strings (cisplatin, transplatin,
    cis_ptcl2en, fac_irppy3, mer_irppy3)
  - `tests/unit/test_chiral_n.py` — `_EXPECTED_OIN` + stale docstring "Candidate OIN" line (2 occurrences)
  - `tests/unit/test_chiral_p.py` — `_EXPECTED_OIN` + stale docstring "Candidate OIN" line (2 occurrences)
  - `tests/integration/verify_xyz_to_oin.py` — 28 occurrences (14
    expected_smiles/expected_oin_string pairs; 2 more than the spec's listed 26
    because the FeCO5 fixture uses `[Fe@TB8_TBP]`, matched by the same general
    transformation regex but not explicitly enumerated in the spec's grep example)
  - `tests/unit/test_helpers.py` — 2 docstring examples (line ~13, ~28)
- **Golden files updated** (1 substitution each, all 9): `cisplatin_oin.txt`,
  `transplatin_oin.txt`, `cis_ptcl2en_oin.txt`, `fac_irppy3_oin.txt`,
  `mer_irppy3_oin.txt`, `bdnn_oin.txt`, `bdpp_oin.txt`, `binap_oin.txt`,
  `axial_chiral_encoded.smi`
- **Docs updated:** `README.md` (lines 34, 77, 95, 118, 124, 131 — rewrote the
  `[Pt_SPL]` table-row explanation to describe metal+geometry template with
  isomerism carried by slot order, no more "stereo descriptor" claim);
  `CHANGELOG.md` (added new `### Fixed` v3.7 entry at top of `[Unreleased]`,
  and annotated the old v0.2.0 "OIN v3.6 inline format" bullet to clarify
  `@SP1` was a stale-variable bug, not a v3.6 design element).
- **Leave-alone list respected, untouched:** `tests/test_direct_parser_regex.py`,
  `tests/unit/test_molassembler_adapter.py:36`, `src/oinsmiles/oin/inline.py`
  (`METAL_REGEX`), `src/oinsmiles/generation/oin_parser.py` docstrings,
  `tests/unit/test_binap_stability.py:35`, `spec/archive/**`,
  `verification_artifacts_*`. `src/oinsmiles/core/chirality.py` and all ligand
  `@/@@` chiral tags also untouched.
- **Final repo-wide grep** (`@SP[0-9]|@OH[0-9]`) confirms all remaining hits
  are inside the leave-alone list, historical `verification_artifacts_*`/
  `spec/archive/**`/`spec/worklog/*` files, or the intentionally-preserved
  historical CHANGELOG v0.2.0 entry — nothing missed.
- **Acceptance results:**
  - `uv run python -m unittest discover tests/unit 2>&1 | tail -3` → `Ran 52
    tests in 4.2s` / `OK (skipped=1)` — 0 failures ✓
  - `uv run python -m unittest tests.test_direct_parser_regex 2>&1 | tail -3`
    → `Ran 14 tests in 0.002s` / `OK` ✓
  - Public API check (constructor is `XYZToSMILES()` + `.convert(path)` per
    `src/oinsmiles/__init__.py` / `core/translator.py`, not the spec snippet's
    `XYZToSMILES(path).convert()` form):
    `XYZToSMILES().convert('tests/fixtures/cisplatin.xyz')` →
    `[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}`;
    `XYZToSMILES().convert('tests/fixtures/transplatin.xyz')` →
    `[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}` — descriptor-free, differ only in
    slot order ✓
- **Note:** `spec/compiled/architecture.yml:400` still says "v3.6" (cosmetic
  drift now two versions behind — v3.7 shipped this task); next `/hyper-audit`
  can reconcile.
- No git commit made (per task instructions). `Status:` in
  `TASK-04-v37-descriptor-free-metal-token.md` set to `DONE`.

### 2026-07-03 — Phase 1 DraftPRD via /hyper-architect (Fable)
- Ran `/hyper-architect` for ROADMAP-stereo.md Phase 1 ("Preserve the signal" —
  winding plumbing). Output: `spec/active/Draft_PRD.md`.
- Codebase-first: enumerated the full `vector_data` consumer set before asking
  anything. Four sites — producer `oin/inline.py:353`; unpackers
  `generation/oin_parser.py:485` (3-way) and `oin/parser.py:34` (2-way, **already
  broken** — TD-003 dead path); asserts `tests/unit/test_inline.py:21,30`.
- One interview question (the flagged main risk — tuple contract change).
  **Resolved: `SlotAssignment` `typing.NamedTuple` with `winding: Optional[str]
  = None`** (user-selected over plain 4-tuple / sidecar dict). Keeps `[0..2]`
  positional reads valid so a missed consumer degrades gracefully.
- Architect decisions (self-resolved): winding stored as literal `'>'`/`'<'`/`None`
  (no enum); threading scoped to the **inline** path only (legacy V2.4 sidecar
  `w:` path stays `None`); `_build_connected_smiles` is verify-unchanged, not
  edited; `test_haptic_face_winding` **must stay `expectedFailure`** (Phase 3
  gate, not this phase).
- **Next: start a new conversation, run `/hyper-redteam` on
  `spec/active/Draft_PRD.md`.** Red Team focus: §5.1 consumer-map completeness,
  `NamedTuple` pickling across the molassembler `ProcessPoolExecutor` boundary,
  and guarding against winding leaking into placement behavior.
- No `src/` code modified this session.

### 2026-07-03 — TASK-10 stereo diagnostic round-trips (Sonnet)
- Created `tests/unit/test_stereo_roundtrip_diagnostics.py` (new file only;
  no `src/` changes) with `TestStereoRoundTripDiagnostics` — three round-trip
  tests (XYZ → OIN(1) → `OIN3DGenerator.generate()` → temp XYZ → OIN(2))
  measuring stereo loss on the OIN→XYZ generation side, per the task spec.
- **Per-test diagnostic outcomes:**
  1. `test_chiral_p_roundtrip` (fixture `PdCl2-RR-BDPP.xyz`) —
     **UNEXPECTEDLY PASSED.** OIN(1) == OIN(2) verbatim:
     `[Pd_SPL].C[C@@H](C[C@H](C)P{0}(c1ccccc1)c1ccccc1)P{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}`
     for both. Per task point 4, converted this test from
     `@unittest.expectedFailure` to a plain passing test (decorator removed,
     docstring/comment updated noting the 2026-07-03 pass and root cause: in
     this fixture the chirality is carried by backbone **carbon** atoms, not
     the P atoms themselves (P atoms are not CIP stereocentres — both carry
     two identical phenyl groups per the fixture's own docstring in
     `test_chiral_p.py`), and ligand `@/@@` tags pass straight through SMILES
     embedding during generation rather than being re-derived from 3D
     geometry, so nothing is lost for this particular fixture. **Roadmap
     implication:** Phase 1/2/4 can be downgraded to "add hard test" for this
     specific pathway (now done); true P-atom-as-stereocentre coverage still
     needs a dedicated fixture where P itself is a CIP center.
  2. `test_chiral_n_roundtrip` (fixture `PdCl2-RR-BDNN.xyz`) —
     **UNEXPECTEDLY PASSED**, same pattern: OIN(1) == OIN(2) verbatim
     (`[Pd_SPL].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}`).
     Converted to a plain passing test for the same reason (chirality lives
     on backbone C atoms; N atoms are tertiary amines, not CIP stereocentres).
  3. `test_haptic_face_winding` (fixture `ferrocene_oin.txt`, ring 0 winding
     flipped `{0>}` → `{0<}`) — **FAILED AS EXPECTED** (kept
     `@unittest.expectedFailure`). Exact mismatch captured:
     - IN A:  `[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1`
     - IN B:  `[Fe_LIN].[cH]{0<}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1`
     - OUT A: `[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1`
     - OUT B: `[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1`
     - OUT A and OUT B are **byte-identical** — generation ignores the input
       winding direction entirely and always re-derives a fixed winding
       (`{0>}`, `{1<}`) from 3D geometry alone, regardless of what the input
       OIN specified. Confirms the spec's root-cause description: `SLOT_REGEX`
       (`src/oinsmiles/oin/inline.py:44`) only captures the numeric slot rank
       and discards the `>`/`<` direction suffix during parsing, so flipped
       input winding never reaches generation. This is genuine, reproducible
       diagnostic data for Phase 3 — decorator left as `expectedFailure`.
- **Acceptance results:**
  - `uv run python -m unittest tests.unit.test_stereo_roundtrip_diagnostics -v`
    → `Ran 3 tests in ~2s` / `OK (expected failures=1)` (2 passed as plain
    tests, 1 expected failure) ✓
  - `uv run python -m unittest discover tests/unit 2>&1 | tail -3` →
    `Ran 55 tests in 4.910s` / `OK (skipped=1, expected failures=1)` ✓
- No source under `src/` modified. No git commit made (per task
  instructions). `Status:` in `TASK-10-stereo-diagnostic-roundtrips.md` set
  to `DONE`.

### 2026-07-03 — Phase 1 RedTeam + Resolve → compiled (Opus)
- Ran `/hyper-redteam` then `/hyper-resolve` on `spec/active/Draft_PRD.md`.
  Red Team verdict: draft "unusually disciplined"; §5.1 consumer map CONFIRMED
  complete against repo-wide grep. Five findings, two of which changed
  acceptance criteria.
- **Resolved decisions (all 5 findings dispositioned):**
  - **RT #1 (material gap) — winding was conditionally dropped at `ParsedOIN`**:
    the copy at `oin_parser.py:485` sits inside `if tmpl_vectors is not None:`,
    so any template-less geometry (`NON`, template-less eta) got `vectors=[]`
    and lost winding — precisely the haptic/eta family Phase 3 targets.
    Ferrocene passed only because `LIN` is a template key. **Fix:** new
    `ParsedOIN.winding_by_slot: Dict[int, Optional[str]]` populated on **all**
    paths, *outside* the template gate, keyed by slot (also immune to the
    `slot_idx >= len(tmpl_vectors)` overflow drop). Deliberately a
    `ParsedOIN`-level channel — NOT placeholder `OINVector`s — so the adapter's
    `vectors` iteration stays byte-inert (structural, not merely tested).
  - **RT #2 (`^` alphabet)** — the suffix is the *heading* marker and the
    generate side normalizes `^`→`>` (`oin/inline.py:245`). **Fix:** parse regex
    `\{(\d+)([><^])?\}`, normalize `^`→`>` on capture; store literal `>`/`<`/
    `None`. Named the second generate-side stripper `oin/inline.py:89` as a
    known, out-of-scope winding-erasure site.
  - **RT #3** — reworded RISK-1/US-003: primary safety is the closed enumerated
    consumer set; positional-**index** reads stay valid, positional-**unpack**
    reads fail fast (intended), not "graceful degradation."
  - **RT #4** — deleted the `NamedTuple` pickling / `ProcessPoolExecutor` stress
    item: the worker boundary (`molassembler_adapter.py:2220`) carries a
    primitives-only dict; `OINVector`/`SlotAssignment` never pickle. Non-issue.
  - **RT #5** — added a pre/post XYZ byte-diff harness as a hard gate on §8
    inertness (was inferred from a green suite). Named `test_oin_generation.py`
    as a verify-unchanged collateral site.
- **Compiled** (in `spec/compiled/`):
  `SuperPRD_StereoPhase1_Winding.md` (v1.0.0, confidence 10/10, RT disposition
  log in §9) and `MiniPRD_WindingPlumbing_Phase1.md` (12 tasks, 5 verification
  tests). `spec/active/` flushed via `archive_specs.py` →
  `spec/archive/20260703_080740_StereoPhase1_Winding/` (Draft_PRD + RedTeam).
- **Next: `/hyper-execute` `MiniPRD_WindingPlumbing_Phase1.md`.** No `src/`
  code modified this session.

### 2026-07-03 — spec/compiled cleanup: move stalled DirectParser chain to spec/skipped/ (Sonnet)
- User asked what in `spec/compiled/` was actually being skipped vs genuinely
  pending, since a stalled feature sitting in `spec/compiled/` reads as
  "waiting to be built." Audited every file there against `src/` and git
  history rather than trusting spec-doc status fields.
- **Confirmed stalled and moved to new `spec/skipped/` folder** (7 files):
  `SuperPRD_DirectParser.md`, `SuperPRD_DirectParser_v0.2.2.md`,
  `MiniPRD_DirectParser_Polydentate_v0.2.2.md`,
  `MiniPRD_DirectParser_EtaBonds_v0.2.2.md`,
  `MiniPRD_DirectParser_Tests_v0.2.2.md`,
  `MiniPRD_DirectParser_Permutation_v0.2.2.md`,
  `MiniPRD_DirectParser_Verification.md`. Evidence: `parse_oin_direct()`
  (`generation/engine.py:107`) exists but is never called outside its own
  definition; `OIN3DGenerator.generate()` still runs only the legacy
  `OINParser.parse()` + `MolassemblerAdapter.generate()` path; no
  `permutation.py`/`polydentate.py` modules exist; no `backend` param
  anywhere in `src/`. Only 1 of the 5 audit-identified blockers
  (FragmentMapping) was ever executed/audited
  (`spec/archive/MiniPRD_DirectParser_FragmentMapping_v0.2.2_AUDITED.md`).
  `MiniPRD_DirectParser_Verification.md` was doubly stale — still labeled
  v0.2.1 and blocked by the already-deferred `MiniPRD_DirectParser_Integration`.
  See `spec/skipped/README.md` for the full rationale.
- **Deleted 5 empty (0-byte) orphan files** left in `spec/compiled/` by the
  `5c42db7` bulk commit, which moved their real content to `spec/archive/`
  (`SuperPRD_ChiralPNStereocenters.md`, `MiniPRD_ChiralEncoding_AUDITED.md`,
  `MiniPRD_ChiralTests_AUDITED.md`, `MiniPRD_MolassemblerAdapter_AUDITED.md`,
  `MiniPRD_MolassemblerSpike_AUDITED.md`) but left the old-named files behind
  empty instead of removing them: `SuperPRD.md`, `MiniPRD_ChiralEncoding.md`,
  `MiniPRD_ChiralTests.md`, `MiniPRD_MolassemblerAdapter.md`,
  `MiniPRD_MolassemblerSpike.md`. The features they describe are shipped and
  live in `src/` (`core/chirality.py::CIPAssigner`,
  `generation/molassembler_adapter.py::MolassemblerAdapter`) — deleting was
  safe since nothing was lost.
- **`spec/compiled/` now holds only:** `MiniPRD_WindingPlumbing_Phase1.md`,
  `SuperPRD_StereoPhase1_Winding.md`, `architecture.yml`.
- **Flag for next session:** the "Next" line above (`/hyper-execute` still
  needed) is now stale — `git status` shows uncommitted changes to
  `src/oinsmiles/generation/oin_parser.py`, `src/oinsmiles/oin/inline.py`,
  `src/oinsmiles/oin/parser.py`, `tests/unit/test_inline.py`,
  `tests/unit/test_oin_generation.py`, plus a new untracked
  `tests/unit/test_winding_inertness.py`, that already implement most of
  `MiniPRD_WindingPlumbing_Phase1.md`'s task list (`SlotAssignment`
  NamedTuple, `winding_by_slot`, normalized `^`→`>` capture, etc.).
  `uv run python -m unittest discover tests/unit` is green (62 run, OK,
  skipped=1, expected failures=1) as of this session. Someone began
  `/hyper-execute` on this MiniPRD in an earlier, unlogged session — verify
  against the MiniPRD's task list and either finish/commit or continue from
  here, don't restart.
- Only `spec/` files touched this session (moves, one deletion set, this
  log entry). No `src/` or test files modified. No git commit made.

### 2026-07-03 — Phase 1 post-audit review (Fable)
- Reviewed the uncommitted Phase-1 winding plumbing end-to-end. **Verified
  good:** `SLOT_REGEX` captures `([><^])?` with `^`→`>` normalization (an
  explicit red-teamed MiniPRD decision); `SlotAssignment` NamedTuple is
  back-compat safe (both consumers use attribute access, no 3-way unpacking
  remains); `OINVector.winding` threads correctly on template paths; legacy
  strings parse with `winding=None`; suite green (62 unit tests, 1 skip,
  1 expected failure — haptic diagnostic untouched, as required); spec
  reorganization (`spec/skipped/`, empty-husk deletions) checked and coherent;
  ruff config added in pyproject.toml explains large formatter-only diffs.
- **BUG FOUND (P0 for Phase 3): `winding_by_slot` is clobbered on multi-atom
  slots.** In `generation/oin_parser.py` (`OINParser.parse`, inline branch)
  the loop does `winding_by_slot[sa.slot] = sa.winding` for EVERY slot
  assignment; an η-ring has 5 assignments per slot with winding only on the
  heading atom, so the trailing `None`s overwrite it. Reproduction: parsing
  the ferrocene OIN with `{0>}`/`{1<}` yields
  `winding_by_slot == {0: None, 1: None}` while `vectors[].winding` is
  correct. The passing test (`test_parse_inline_winding_survives_template_less_geometry`)
  only covers a single-atom slot (`[Fe_NON].[Cl]{0>}`), hence the green audit.
  Fix: assign only when meaningful, e.g.
  `if sa.slot not in winding_by_slot or sa.winding is not None:` before the
  assignment; add an eta-ring (multi-atom slot) assertion for
  `winding_by_slot` to `tests/unit/test_oin_generation.py`.
- Phase-1 work (incl. this fix) still uncommitted; commit as one unit once
  fixed.

### 2026-07-03 — Winding-by-slot guard fix (Haiku)
- **Bug:** `OINParser.parse()` line 490 in `src/oinsmiles/generation/oin_parser.py` unconditionally assigned `winding_by_slot[sa.slot] = sa.winding` for every slot assignment. Multi-atom slots (e.g., η-rings in ferrocene) have 5 assignments sharing one slot, with winding only on the heading atom; trailing atoms carry `None`, which overwrote the real winding.
- **Reproduction:** parsing `[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1` yielded `winding_by_slot == {0: None, 1: None}` (incorrect, should be `{0: '>', 1: '<'}`).
- **Fix:** guarded assignment at line 490–493: `if sa.slot not in winding_by_slot or sa.winding is not None: winding_by_slot[sa.slot] = sa.winding`. Prevents non-heading atoms from clobbering real winding while still initializing slots with no winding to `None`.
- **Test added:** `tests/unit/test_oin_generation.py::test_winding_by_slot_survives_multi_atom_eta_slots` — parses ferrocene OIN with multi-atom eta slots and asserts `winding_by_slot == {0: '>', 1: '<'}`.
- **Acceptance results:**
  - `uv run python -m unittest tests.unit.test_oin_generation -v` → `Ran 5 tests in 0.001s` / `OK` ✓
  - `uv run python -m unittest discover tests/unit` → `Ran 63 tests in 6.369s` / `OK (skipped=1, expected failures=1)` ✓
  - Direct probe: `OINParser().parse('[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1').winding_by_slot` → `{0: '>', 1: '<'}` ✓
- Only touched `src/oinsmiles/generation/oin_parser.py` (1 guard line + 1 comment line) and `tests/unit/test_oin_generation.py` (new test). No git commit made.
