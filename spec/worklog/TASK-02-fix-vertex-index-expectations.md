# TASK-02: Fix vertex_indices expectations to 1-based fragment ranks

Status: DONE
Depends on: none
Suggested model: Haiku

## Goal

7 tests in `tests/test_direct_parser_regex.py` fail because they expect
`vertex_indices` to be 0-based, but the implementation correctly returns
**1-based fragment ranks**. Fix the test expectations, NOT the implementation.

## Context (no prior repo knowledge needed)

`_extract_oin_constraints()` in `src/oinsmiles/generation/oin_parser.py`
parses OIN strings like `[Pt@SP1_SQP].[Cl]{0}.[Cl]{1}` and returns constraint
dicts whose `vertex_indices` values are **fragment ranks**: the metal is
always fragment rank 0, and ligand fragments are enumerated starting at 1
(see `src/oinsmiles/generation/engine.py:158` — `updated_frag_to_atom[0] = [0]`
for the metal, ligands `enumerate(..., start=1)`). Downstream,
`construct_molassembler_mol` (`src/oinsmiles/generation/oin_parser.py:153-167`)
indexes `frag_to_atom[frag_rank]`; a 0-based ligand value would resolve to the
metal atom itself and trip the bond-to-self guard at
`src/oinsmiles/generation/oin_parser.py:180-184`. So 1-based is the correct,
load-bearing convention (decision D-2 in `spec/worklog/NOTES.md`).

The list POSITION corresponds to the slot; the VALUE is the fragment rank.

## Files to touch (only this one)

- `tests/test_direct_parser_regex.py`

## Steps

Update the expected lists at these lines (line numbers as of 2026-07-02):

| Line | Test | Current expectation | Change to |
|---|---|---|---|
| 20 | test_platinum_square_planar_basic | `[0, 1, 2, 3]` | `[1, 2, 3, 4]` |
| 42 | test_no_shape_code | `[0, 1]` | `[1, 2]` |
| 58 | test_no_chiral_tag | `[0, 1]` | `[1, 2]` |
| 71 | test_iron_linear_with_heading_markers | `[0, 1]` | `[1, 2]` |
| 112 | test_heading_marker_clockwise | `[0, 1]` | `[1, 2]` |
| 123 | test_heading_marker_counterclockwise | `[0, 1]` | `[1, 2]` |
| 132 | test_multiple_vertex_indices (`expected_indices = ...`) | `[0, 1, 2, 3, 4, 5]` | `[1, 2, 3, 4, 5, 6]` |

Also add one comment near the first fixed assertion documenting the
convention, e.g.:
`# vertex_indices values are fragment ranks (metal = rank 0, ligands start at 1); list position = slot`

## Acceptance (exact commands)

```
uv run python -m unittest tests.test_direct_parser_regex -v 2>&1 | tail -5
```
Expected: **Ran 14 tests ... OK** (was 7 failures).

## Constraints / DO NOT

- Do NOT modify `src/oinsmiles/generation/oin_parser.py` or any other source
  file. If a test still fails after the expectation change, STOP and report —
  do not adjust the implementation.
- Do NOT change the literal `@SP1` inputs in the test strings — they are
  deliberate back-compat coverage for legacy-format parsing (decision D-3).

## Out of scope

- Everything else in the suite (TASK-01/03/04).

## On completion

Set `Status: DONE` above and append a dated Log entry to
`spec/worklog/NOTES.md`.
