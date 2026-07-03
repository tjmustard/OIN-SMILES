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

**Current status (2026-07-03, after TASK-01/02/03): ALL GREEN.**
- `uv run python -m unittest discover tests` → 55 run, OK
- `uv run python -m unittest discover tests/unit` → 52 run, OK (skipped=1)
- Note: root discovery does not recurse into `tests/unit` (no `__init__.py`
  there — see AGENTS.md); run both commands to cover everything.
- Changes uncommitted as of this entry. Next up: TASK-04 (v3.7 revert).

| Task | Status | Model tier | Depends on |
|---|---|---|---|
| TASK-01 delete stale root tests | DONE | Haiku | — |
| TASK-02 vertex-index expectations | DONE | Haiku | — |
| TASK-03 GeneratedStructure assert | DONE | Haiku | — |
| TASK-04 v3.7 descriptor-free metal token | TODO | Sonnet | TASK-01 |
| TASK-10 stereo diagnostic round-trips | TODO | Sonnet | TASK-04 |
| Stereo Phases 1–4 | ROADMAP (needs MiniPRDs) | Sonnet/Opus | TASK-10 |

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
