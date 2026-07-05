# TASK-41: E3 — raise a real error instead of bare `None` (WS-1)

Status: DONE (2026-07-04) — bare `None`→descriptive `ValueError` at the
`get_lig_mol` guard; TiCat3/4 now report `get_lig_mol failed for ligand fragment
#1 (SMILES: ...)`. Fixture `tests/fixtures/ticat3_generated_broken.xyz` +
`tests/unit/test_xyz2mol_errors.py` pin it (deterministic, reproduces). Suites:
`discover tests` 55 OK, `discover tests/unit` 125 OK/3 skip/0 xfail,
`verify_xyz_to_oin.py` 25/25. See NOTES.md Log.
Depends on: TASK-40 (uses `--only` iteration + the failure artifact)
Suggested model: Haiku (source fix is trivial; fixture-capture steps are exact)
Effort: ~0.25 session · Risk: none (touches an already-crashing path only)

Part of the η-ligand round-trip recovery effort — see
`spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` (§2 root cause E3; §3 WS-1).

## Goal

`get_tmc_mol` returns a bare `None` when ligand perception fails. Its sole
convert-path caller unpacks a 2-tuple, so the failure surfaces as the opaque
`cannot unpack non-iterable NoneType object` (TiCat3/4). Replace the bare
`return None` with a descriptive `ValueError` naming the offending fragment. This
flips nothing green — it makes the existing crash diagnosable.

## Root cause (exact)

`src/oinsmiles/utils/xyz2mol.py`, function `get_tmc_mol` (def ~line 472):

```python
lig_mol, lig_charge = get_lig_mol(m, lig_charge, lig_coordinating_atoms)
if not lig_mol:
    return None            # <-- bare None; get_tmc_mol otherwise returns a 2-tuple
```

Caller `src/oinsmiles/core/translator.py:27-30`:

```python
try:
    tmc_mol, xyz_coords = get_tmc_mol(path, charge, with_stereo=False)  # unpack crash on None
except Exception as e:
    raise ValueError(f"xyz2mol failed: {e}")
```

## Audit (confirmed 2026-07-04 — record and respect)

The **only** bare-`None` return in `get_tmc_mol` is this one (~line 543). The two
`return None, charge` at `xyz2mol.py:372` and `:387` are inside a **different**
function, `get_lig_mol`, whose contract is a 2-tuple `(mol, charge)` that the
caller already handles at the `if not lig_mol` guard. They are NOT
unpacking bugs — **leave them alone.**

## The fix — one file, one return site: `src/oinsmiles/utils/xyz2mol.py`

Replace the `return None` in the `if not lig_mol:` guard with:

```python
if not lig_mol:
    try:
        frag_smiles = Chem.MolToSmiles(m)
    except Exception:
        frag_smiles = f"<{m.GetNumAtoms()} atoms>"
    raise ValueError(
        f"get_lig_mol failed for ligand fragment #{i} "
        f"(SMILES: {frag_smiles!r}); cannot build TMC mol"
    )
```

`translator.py` already wraps this, so the message surfaces as
`xyz2mol failed: get_lig_mol failed for ligand fragment #N (SMILES: '...'); ...`.
No caller change is needed (grep-confirmed sole convert-path caller).

## New unit test (deterministic fixture)

The TiCat3 *generated* XYZ is written by the harness at step 2 **before** the
step-3 re-encode crash (ETKDG uses a fixed seed → reproducible). Capture it as a
committed fixture and pin the error:

1. `uv run python tests/integration/verify_roundtrip.py --only TiCat3 --output-dir /tmp/rt41`
2. Copy the generated structure to a fixture:
   `cp /tmp/rt41/Ex21_TiCat3_generated.xyz tests/fixtures/ticat3_generated_broken.xyz`
   (confirm the `Ex21_..._generated.xyz` file exists first; adjust the index if
   the example ordering differs).
3. Add `tests/unit/test_xyz2mol_errors.py`:
   ```python
   import unittest
   from pathlib import Path

   from oinsmiles.utils.xyz2mol import get_tmc_mol


   class TestXyz2MolErrors(unittest.TestCase):
       def test_get_tmc_mol_raises_valueerror_on_unbuildable_ligand(self):
           fx = Path(__file__).parent.parent / "fixtures" / "ticat3_generated_broken.xyz"
           with self.assertRaises(ValueError) as ctx:
               get_tmc_mol(fx, 0, with_stereo=False)
           self.assertIn("get_lig_mol failed", str(ctx.exception))
   ```
   If the fixture proves non-reproducible (get_tmc_mol does NOT fail on it),
   `@unittest.skipUnless`-guard the test and report in the Log — do NOT force it.

## Acceptance (exact commands + expected)

1. Descriptive message on the crashing complexes:
   ```
   uv run python tests/integration/verify_roundtrip.py --only TiCat3 --output-dir /tmp/rt41
   cat /tmp/rt41/Ex*_TiCat3_error.txt | head -3
   ```
   Expect: message contains `get_lig_mol failed for ligand fragment #` (no longer
   `cannot unpack non-iterable NoneType`). TiCat3/4 still FAIL overall.
2. New unit test green:
   ```
   uv run python -m unittest tests.unit.test_xyz2mol_errors -v 2>&1 | tail -3
   ```
3. Full suites green, no regression:
   ```
   uv run python -m unittest discover tests 2>&1 | tail -3          # 55 OK
   uv run python -m unittest discover tests/unit 2>&1 | tail -3     # 125 OK, skipped=3, xfail=0
   uv run python tests/integration/verify_xyz_to_oin.py 2>&1 | tail -2  # 25 Passed
   ```

## Constraints / DO NOT

- Touch only `src/oinsmiles/utils/xyz2mol.py` (one return site),
  `tests/unit/test_xyz2mol_errors.py` (new), and the new fixture.
- Do NOT touch the `get_lig_mol` 2-tuple returns (`:372`, `:387`) or any
  currently-passing path.
- Run `uv run ruff format` on the edited source + test files. Re-read after
  autofix (ruff-docstring-truncation gotcha).

## On completion

Set `Status: DONE`, append a dated Log entry to `spec/worklog/NOTES.md`, update
the handoff status table (WS-1 landed; the unit count is now 125). Do NOT commit
unless asked — leave staged for review (scoped `git add`).
