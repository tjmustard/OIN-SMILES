# TASK-04: OIN v3.7 — descriptor-free metal token (fix stale `is_metal` bug)

Status: DONE
Depends on: TASK-01 (the deleted root tests assert the old style)
Suggested model: Sonnet

## Goal

Restore the intended OIN string style: metal token `[Pt_SPL]` instead of
`[Pt@SP1_SPL]`. The `@SP1`/`@OH10` segment is NOT a designed feature — it is
RDKit's non-tetrahedral SMILES extension leaking through a stale-variable bug,
and it is redundant (cis/trans and fac/mer are fully encoded by slot ordering)
and atom-order dependent (README admits this). This change is labeled **OIN
v3.7** in docs/CHANGELOG. Parsers already tolerate and discard `@desc`, so old
strings keep working.

## Context (no prior repo knowledge needed)

Emission chain of the bug (verified 2026-07-02, line numbers current then):

1. `src/oinsmiles/core/chirality.py:98` — `Chem.AssignAtomChiralTagsFromStructure(mol)`
   sets `CHI_SQUAREPLANAR`/`CHI_OCTAHEDRAL` on the metal atom from 3D (this
   call is required for P/N tetrahedral tags — the metal tag is a side effect).
2. `src/oinsmiles/utils/xyz2mol.py:868` — first per-fragment loop sets local
   `is_metal = (metal_idx in indices)` and stores it in
   `fragments_data[i]['is_metal']` (line 901).
3. **THE BUG**: the second per-fragment loop (`for i, item in
   enumerate(fragments_data):` at line 910) never reassigns `is_metal`. At
   line 982 (`if not is_metal:`) it holds the stale value from the LAST
   fragment of loop 1 (a ligand → False), so the METAL fragment takes the
   ligand branch (lines 983–986) where `Chem.MolToSmiles(...,
   isomericSmiles=True)` serializes the metal chiral tag as `[Pt@SP1]`.
   The intended metal branch at line 988
   (`sanitized_smiles = f"[{...GetSymbol()}]"`) is dead code in practice.
4. `src/oinsmiles/oin/inline.py:105` — `generate_inline_string` appends
   `_GEO` → `[Pt@SP1_SPL]`.

### The fix (one line)

In `src/oinsmiles/utils/xyz2mol.py`, at the top of the second loop body
(immediately after line 910's `for i, item in enumerate(fragments_data):`,
next to the existing `indices = item['indices']`), add:

```python
is_metal = item['is_metal']
```

The metal fragment then flows through line 988 → `[Pt]` → `[Pt_SPL]`.

Note this also fixes a latent inverse hazard: if the metal were ever the LAST
fragment, the old code emitted `[El]` for EVERY fragment.

### String transformation rule for expectations/goldens

Remove only the metal token's `@DESC` segment (uppercase letters+digits
between element and `_GEO`):
regex `\[([A-Z][a-z]?)@[A-Z]+[0-9]+_` → `[\1_`
Examples: `[Pt@SP1_SPL]`→`[Pt_SPL]`, `[Pt@SP2_SPL]`→`[Pt_SPL]`,
`[Ir@OH10_OCT]`→`[Ir_OCT]`, `[Pd@SP1_SPL]`→`[Pd_SPL]`,
`[Fe@OH19_OCT]`→`[Fe_OCT]`, `[V@OH21_SPY]`→`[V_SPY]`.
**NEVER touch ligand chiral tags like `[C@@H]`/`[C@H]`** — they have no
`_GEO` suffix and the regex above cannot match them (they lack the digit+`_`).
Note the transformation is lossy across isomer pairs (`@SP1` and `@SP2` both
→ `[Pt_SPL]`): cis vs trans stays distinguished by the slot ordering that
follows, e.g. cis `.[Cl]{0}.[Cl]{1}.N{2}.N{3}` vs trans
`.[Cl]{0}.N{1}.[Cl]{2}.N{3}`.

## Files to touch

**Source (the fix):**
- `src/oinsmiles/utils/xyz2mol.py` — the one line above.

**Test expectations (apply transformation rule):**
- `tests/unit/test_regression_stability.py` lines 31, 37, 43, 55, 61
- `tests/unit/test_chiral_n.py` line 33 (`_EXPECTED_OIN`, `[Pd@SP1_SPL]`)
- `tests/unit/test_chiral_p.py` line 34 (`_EXPECTED_OIN`, `[Pd@SP1_SPL]`)
- `tests/integration/verify_xyz_to_oin.py` — 26 occurrences (13
  expected_smiles/expected_oin_string pairs) at lines 263-264, 279-280,
  295-296, 327-328, 342-343, 357-358, 372-373, 387-388, 433-434, 510-511 and
  any others matching the pattern (grep to be exhaustive; includes `@SP3`,
  `@OH19`, `@OH21`)

**Golden files (apply transformation rule to each single-line file):**
- `tests/candidate_outputs/cisplatin_oin.txt`, `transplatin_oin.txt`,
  `cis_ptcl2en_oin.txt`, `fac_irppy3_oin.txt`, `mer_irppy3_oin.txt`,
  `bdnn_oin.txt`, `bdpp_oin.txt`, `binap_oin.txt`, `axial_chiral_encoded.smi`

**Docs:**
- `tests/unit/test_helpers.py` lines 13, 28 (docstring examples)
- `README.md` lines 34, 77, 95, 118, 124, 131 — line 131's table row currently
  explains `SP1` as "stereo descriptor (atom-order dependent)"; rewrite the row
  for `[Pt_SPL]`: metal + geometry template, isomerism carried by slot order.
- `CHANGELOG.md` — add a new **v3.7** entry at top: "OIN v3.7: descriptor-free
  metal token (`[Pt_SPL]`, was `[Pt@SP1_SPL]`). The `@desc` was an RDKit
  non-tetrahedral stereo leak via a stale `is_metal` variable in xyz2mol.py;
  isomer information was and remains fully encoded by slot ordering. Parsers
  continue to accept legacy `@desc` strings." (Also fix line 38's claim that
  the `@` form is "v3.6" if editing nearby.)

**Leave alone (deliberate legacy inputs / tolerant parsers):**
- `tests/test_direct_parser_regex.py` — feeds literal `@SP1` INPUT strings;
  this is the back-compat parsing coverage. Do not change.
- `tests/unit/test_molassembler_adapter.py:36` — `original_oin` input; parser
  tolerance coverage. Do not change.
- `src/oinsmiles/oin/inline.py` (`METAL_REGEX` :41 and comments) and
  `src/oinsmiles/generation/oin_parser.py` (docstrings :220, :230) — the
  tolerant parsing must remain; docstring wording update optional.
- `tests/unit/test_binap_stability.py:35` — comment only; optional cosmetic.
- `spec/archive/**`, `verification_artifacts_*` — historical; do not touch.

## Steps

1. Apply the one-line source fix.
2. Run `uv run python -m unittest discover tests/unit` — expect the
   regression/chiral tests to now FAIL against old expectations (proves the
   fix emits descriptor-free strings).
3. Update all expectations/goldens/docs per the lists above (grep
   `@SP[0-9]\|@OH[0-9]` repo-wide afterward to confirm nothing missed outside
   the leave-alone list).
4. Run acceptance.

## Acceptance (exact commands)

```
uv run python -m unittest discover tests/unit 2>&1 | tail -3
```
Expected: 52 run, **0 failures**, 1 skip (assumes TASK-03 also done; if not,
the only failure is TASK-03's known one).

```
uv run python -m unittest tests.test_direct_parser_regex 2>&1 | tail -3
```
Expected: OK (legacy `@desc` inputs still parse — back-compat intact).

```
uv run python -c "
from oinsmiles import XYZToSMILES
print(XYZToSMILES('tests/fixtures/cisplatin.xyz').convert())
print(XYZToSMILES('tests/fixtures/transplatin.xyz').convert())
"
```
Expected: `[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}` and
`[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}` (check the exact public API signature in
`src/oinsmiles/__init__.py` / README if this snippet's constructor form is
wrong; the important part is cis and trans strings differ only in slot order,
with no `@`).

## Constraints / DO NOT

- Do NOT modify `src/oinsmiles/core/chirality.py` — `CIPAssigner`'s
  tetrahedral P/N tags are load-bearing; the metal-side effect is now harmless
  because the metal branch emits a bare `[El]`.
- Do NOT change `METAL_REGEX` or any parser to reject `@desc`.
- Do NOT touch ligand `@/@@` tags anywhere.
- If after the fix any string differs from expectation by MORE than the
  removed `@desc` (e.g. slot order changed), STOP and report — that would
  contradict the verified analysis.

## Out of scope

- Winding markers, P/N stereo on generation, haptic faces (ROADMAP-stereo.md).

## On completion

Set `Status: DONE` above, append a dated Log entry to `spec/worklog/NOTES.md`,
and note that `spec/compiled/architecture.yml:400` still says "v3.6" (cosmetic,
next `/hyper-audit` can reconcile).
