# TASK-01: Delete stale root-level duplicate test files

Status: DONE
Depends on: none
Suggested model: Haiku

## Goal

Remove five outdated test files at `tests/` root (plus their orphaned helper)
that are older duplicates of the maintained copies in `tests/unit/`. They use
package-relative imports (`from .test_helpers import ...`) that break under
`unittest discover` from the repo root, producing 5 import errors, and they
assert the legacy `[Pt@SP1_SPL]` string style that TASK-04 removes.

## Context (no prior repo knowledge needed)

- The canonical unit tests live in `tests/unit/` and use absolute imports
  (`from tests.unit.test_helpers import ...`). They are supersets of the root
  copies (e.g. `tests/unit/test_axial_chiral.py` has 3 tests vs 2 in the root
  copy) and all pass.
- The root copies pass when run directly as a package module but ERROR under
  `discover` because they load as top-level modules where relative imports
  have no parent package.
- `tests/test_helpers.py` (root) is imported ONLY by the five files being
  deleted (verified 2026-07-02 via grep), so it becomes orphaned. Re-verify at
  implementation time: `grep -rn "test_helpers" tests/ --include="*.py" | grep -v tests/unit`
  must show no importers other than the deleted files.

## Files to touch (delete only; nothing else)

- `tests/test_axial_chiral.py`
- `tests/test_binap_stability.py`
- `tests/test_chiral_n.py`
- `tests/test_chiral_p.py`
- `tests/test_regression_stability.py`
- `tests/test_helpers.py` (only after re-verifying it is orphaned)

## Steps

1. Run the grep above to confirm `tests/test_helpers.py` has no other importers.
2. `git rm` the six files.
3. Run acceptance.

## Acceptance (exact commands)

```
uv run python -m unittest discover tests 2>&1 | tail -5
```
Expected: **0 errors** (the 5 `ImportError`/`attempted relative import` errors
are gone). Failures unrelated to imports may remain (they belong to TASK-02/03,
7 failures + 1 unit failure as of baseline) — that is OK for this task.

```
uv run python -m unittest discover tests/unit 2>&1 | tail -3
```
Expected: unchanged from baseline (52 run, 1 failure, 1 skip) — proves the
maintained copies were untouched.

## Constraints / DO NOT

- Do NOT delete or modify anything under `tests/unit/`, `tests/integration/`,
  or `tests/candidate_outputs/`.
- Do NOT delete `tests/test_direct_parser_regex.py`,
  `tests/test_fragment_mapping.py`, `tests/test_oin_parser_ast.py`,
  `tests/test_molassembler_instantiation.py`, or `tests/spike_molassembler.py`
  — those are current, not duplicates.
- Do NOT "fix" the relative imports instead of deleting; deletion is the
  decided action (see spec/worklog/NOTES.md D-1/D-3 context).

## Out of scope

- The 7 `vertex_indices` failures (TASK-02), the adapter assert (TASK-03),
  the `@desc` string style (TASK-04).

## On completion

Set `Status: DONE` above and append a dated Log entry to
`spec/worklog/NOTES.md` (files deleted, acceptance output summary).
