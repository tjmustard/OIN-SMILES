# TASK-03: Update adapter test assertion for GeneratedStructure return type

Status: DONE
Depends on: none
Suggested model: Haiku

## Goal

One test predates the v0.2.0 API change where
`MolassemblerAdapter.generate()` stopped returning a raw XYZ string and now
returns a `GeneratedStructure(xyz: str, mol: Optional[Chem.Mol])` dataclass.
Update the assertion to compare against `result.xyz`.

## Context (no prior repo knowledge needed)

- Failing test: `tests/unit/test_molassembler_adapter.py` →
  `test_successful_generate_returns_xyz_block` (line 115; assertion at 133:
  `self.assertEqual(result, expected_xyz)`).
- Current failure: `GeneratedStructure(xyz='9\n\nPt 0.0 0.0 0.0\n...', mol=None) != '9\n\nPt 0.0 0.0 0.0\n...'`.
- The API change is intentional (see CHANGELOG / `GeneratedStructure` in
  `src/oinsmiles/generation/molassembler_adapter.py`); `.xyz` always holds the
  XYZ block string.

## Files to touch (only this one)

- `tests/unit/test_molassembler_adapter.py`

## Steps

1. Line 133: change `self.assertEqual(result, expected_xyz)` to
   `self.assertEqual(result.xyz, expected_xyz)`.
2. Optionally update the test docstring (line 116) to say the adapter returns
   a `GeneratedStructure` whose `.xyz` is the XYZ block.

## Acceptance (exact commands)

```
uv run python -m unittest tests.unit.test_molassembler_adapter -v 2>&1 | tail -5
```
Expected: all tests **OK** (was 1 failure).

```
uv run python -m unittest discover tests/unit 2>&1 | tail -3
```
Expected: 52 run, 0 failures, 1 skip (the skip is
`test_axial_chiral` — intentional, do not unskip).

## Constraints / DO NOT

- Do NOT modify `src/oinsmiles/generation/molassembler_adapter.py` — the
  return type is correct; only the test is stale.
- Do NOT unskip or edit `tests/unit/test_axial_chiral.py`.

## Out of scope

- TASK-01/02/04.

## On completion

Set `Status: DONE` above and append a dated Log entry to
`spec/worklog/NOTES.md`.
