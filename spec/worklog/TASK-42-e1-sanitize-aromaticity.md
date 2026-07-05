# TASK-42: E1 — sanitize / aromaticity-normalize before serialize (WS-2)

Status: DONE (2026-07-04) — **fix differs from the handoff's proposal; probe
corrected it.** Probe found the handoff's staged sanitize does NOT work: (a)
full `SanitizeMol` raises `KekulizeException` on a charge-less Cp anion, and (b)
`SANITIZE_ALL ^ KEKULIZE` + `SetAromaticity` re-flags atoms but leaves bond
TYPES `SINGLE`, so `[cH]-[cH]` dashes survive. Actual fix = restore
`AROMATIC` bond type on `SINGLE` bonds between aromatic atoms, guarded by
`bond.IsInRing()` (protects biaryl bonds: Ir(ppy)3 phenyl-pyridine, BINAP
binaphthyl). In `generate_robust_smiles` only; `xyz2mol.py:918` NOT touched
(probe showed the in-function fix suffices). Diff gate: 25 pass-1 strings
byte-identical. Ferrocene round-trip PASS (string + RMSD 0.977); full round-trip
19→20/25, zero regressions. TiCp2Me2 string PASS, RMSD 1.675 still FAIL (Phase
1-2 + WS-7). New `tests/unit/test_oin_sanitizer_aromaticity.py` (2 tests) pins
the fix + biaryl guard. Suites: `discover tests` 55 OK, `discover tests/unit`
127 OK/3 skip/0 xfail, `verify_xyz_to_oin.py` 25/25. See NOTES.md Log.
Depends on: TASK-40, TASK-41
Suggested model: **Sonnet** (needs RDKit judgment) — the Phase-0 probe is
load-bearing; do not skip it.
Effort: ~1 session · Risk: charge-less Cp fragments hit kekulize failure →
fallback (b) is load-bearing.

Part of the η-ligand round-trip recovery effort — see
`spec/worklog/ROUNDTRIP-eta-recovery-handoff.md` (§2 root cause E1; §3 WS-2;
§2 feasibility finding #4 — the staged fallback). Read that + `NOTES.md` first.

## Goal

`OINSanitizer.generate_robust_smiles` (`src/oinsmiles/utils/oin_aligner.py:24-75`)
serializes whatever atom flags + bond types are present — it never runs
`SanitizeMol`/`Kekulize`/aromaticity perception. The generator's mol has
aromatic-flagged Cp atoms joined by **SINGLE** bonds, so RDKit emits `[cH]-[cH]`.
This is the primary cause of the Ferrocene/TiCp2Me2 string mismatch on the
step-3 re-encode of the generator's bonded mol. Fix: perceive aromaticity
(with a staged kekulize fallback) before serializing.

## MANDATORY probe first (do not edit source until this is done)

Follow the TASK-30 diagnostic pattern. Write a throwaway script in the
scratchpad (NOT the repo). Measure, record findings in the Log, THEN implement.

1. Reproduce the exact fragment state. Run
   `verify_roundtrip.py --only Ferrocene --output-dir /tmp/probe42`, then in a
   script convert the Ferrocene XYZ → OIN via `OIN3DGenerator().generate(...)`
   and reach into the Cp fragment mol handed to `generate_robust_smiles`
   (or reconstruct it from `gen_result.mol`). Dump, per Cp fragment:
   per-atom `GetIsAromatic()` and per-bond `GetBondType()`. Confirm the
   diagnosis: aromatic atoms joined by SINGLE bonds.
2. On an `RWMol` copy of that fragment, test each staged path and record which
   the 25 complexes' ligand fragments actually need:
   - (a) `Chem.SanitizeMol(rw_mol)` — succeeds, or raises a kekulize error?
   - (b) on kekulize failure:
     `Chem.SanitizeMol(rw_mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)`
     then `Chem.SetAromaticity(rw_mol)`.
   - (c) no-op (current behavior).
   Charge-less Cp anions are expected to require (b).
3. **Decide the second question:** does the fragment-copy bond loop in
   `get_oin_string` (`xyz2mol.py:910-918`, which copies `bond.GetBondType()`
   verbatim) ALSO need bond-type normalization, or is sanitizing inside
   `generate_robust_smiles` sufficient? Measure where the aromatic-single-bond
   state originates (fragment copy vs upstream generator mol).

## The fix

In `src/oinsmiles/utils/oin_aligner.py`, `generate_robust_smiles` (~lines 24-75),
after `rw_mol = Chem.RWMol(ligand_mol)` and BEFORE the H-locking loop (~line 39),
add the staged normalization:

```python
# Perceive aromaticity/valence before serializing so aromatic-flagged atoms
# joined by SINGLE bonds (from the generator's de-aromatized mol) do not emit
# explicit single bonds like [cH]-[cH]. Staged fallback protects charge-less
# Cp-anion fragments that fail kekulization.
try:
    Chem.SanitizeMol(rw_mol)
except Exception:
    try:
        Chem.SanitizeMol(rw_mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        Chem.SetAromaticity(rw_mol)
    except Exception:
        pass  # (c) preserve current no-op serialize behavior
```

Add the `xyz2mol.py:918` bond-type normalization ONLY if the probe (step 3) says
sanitizing inside `generate_robust_smiles` is not sufficient; keep it guarded.

## Strict acceptance (exact commands + expected)

1. **Diff gate — load-bearing.** Pass-1 OIN strings for all 25 complexes
   byte-identical before/after. This protects the re-pinned haptic goldens
   (`test_haptic_face_golden_match` in `tests/unit/test_oin_generation.py`) and
   `tests/unit/test_winding_inertness.py`.
   ```
   # BEFORE editing: uv run python tests/integration/verify_xyz_to_oin.py --output-dir /tmp/e1_before
   # AFTER editing:  uv run python tests/integration/verify_xyz_to_oin.py --output-dir /tmp/e1_after
   diff <(grep 'OIN:' /tmp/e1_before/integration_log.txt) <(grep 'OIN:' /tmp/e1_after/integration_log.txt)
   ```
   Expect: empty diff (still 25/25, emitted strings unchanged).
   **Any diff = regression → STOP, escalate to a MiniPRD. Do not proceed.**
2. Ferrocene round-trip → **PASS** (string + RMSD; its RMSD 0.977 is already
   < 1.0, so E1 flips it fully green):
   ```
   uv run python tests/integration/verify_roundtrip.py --only Ferrocene 2>&1 | grep -E "PASS|FAIL"
   ```
3. TiCp2Me2 string check passes; RMSD stays FAIL (that is Phase 1–2 + WS-7, not
   this task):
   ```
   uv run python tests/integration/verify_roundtrip.py --only TiCp2Me2 2>&1 | grep -E "OIN Stability|RMSD"
   ```
   Expect: `[PASS] OIN Stability`, `[FAIL] Geometry: High RMSD`.
4. New unit test `tests/unit/test_oin_sanitizer_aromaticity.py`: feed
   `generate_robust_smiles` a mol with aromatic-atom + SINGLE-bond Cp state,
   assert no explicit `-` between aromatic ring atoms / consistent aromaticity;
   **pin fallback (b)** with a charge-less Cp-anion fragment (assert it
   serializes aromatic and does not crash or fall to the no-op).
5. Full suites green:
   ```
   uv run python -m unittest discover tests 2>&1 | tail -3          # 55 OK
   uv run python -m unittest discover tests/unit 2>&1 | tail -3     # 126 OK, skipped=3, xfail=0
   uv run python tests/integration/verify_xyz_to_oin.py 2>&1 | tail -2  # 25 Passed
   ```

## Constraints / DO NOT

- Touch `src/oinsmiles/utils/oin_aligner.py` (the sanitize staging), possibly
  `src/oinsmiles/utils/xyz2mol.py:918` (ONLY if the probe requires it, guarded),
  and a new unit test file.
- Do NOT change the slot-mapping logic, the H-locking loop's intent, or the
  `MolToSmiles(..., isomericSmiles=True, canonical=True)` flags.
- Run `uv run ruff format` on edited files; re-read after autofix
  (ruff-docstring-truncation gotcha bit us at commit 033f1c5).

## On completion

Set `Status: DONE`, append a dated Log entry to `spec/worklog/NOTES.md` (probe
findings — which staged path each fragment needed and the step-3 decision — the
fix, the diff-gate result, acceptance results), update the handoff status table
(WS-2 landed; round-trip 19→20/25, Ferrocene green; unit count 126). Do NOT
commit unless asked — leave staged for review (scoped `git add`).
