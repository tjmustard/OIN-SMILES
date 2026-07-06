# Handoff — BDNN Pd square-plane stability (ensemble size)

**Difficulty:** LOW. Suitable for Haiku/Sonnet to run + verify; escalate to Opus/Sonnet
only if the ensemble flag alone proves unreliable.
**Branch:** `feature/metallogen-3d-generator` (do NOT push; commit with `git commit --no-verify`).
**Status when written (2026-07-05):** already empirically fixed by `--ensemble-size 5`
(see `verification_artifacts_BDNN_20260705_185631`). This task is to make it robust and
land it, not to discover the fix.

## Context / problem

Round-trip: `XYZToSMILES(input) → OIN → OIN3DGenerator(engine="metallogen") → XYZ → OIN`,
compared to the input OIN. In the `--ensemble-size 1` full MACE run
(`verification_artifacts_20260705_162501`, roundtrip 19/25), **PdCl2-RR-BDNN failed**:

```
Exp: [Pd_SPL].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}
Got: [Pd_TPY].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}
High RMSD: 1.5681
```

- The **stereo string is correct** (`C[C@@H](C[C@H]...)` in both) — the stereo-carry fix
  (`56a543a`) works for BDNN.
- The failure is purely **geometry**: the generated Pd center distorted from square-planar
  (`SPL`) toward trigonal-pyramidal (`TPY`), giving RMSD 1.57 and a geo-code mismatch.
- It is **stochastic**: BDNN's *dedicated* single-complex MACE run passed at RMSD 0.117
  (clean `SPL`). N-donors (this diamine) are floppier than P-donors, so the square-planar
  Pd is a shallow minimum the embed sometimes misses. BDPP (the P analogue) passed.

## Root cause & why ensemble helps

`generate_3d_structures` (`src/oinsmiles/generator3d/__init__.py`) builds a pool of
conformers over `scales × options`, FF-cleans each, then (if an optimizer is set)
energy-ranks and returns the lowest. With `ensemble_size=1` only one survivor is kept, so a
single bad (TPY-distorted) embed is returned. With `ensemble_size>1` the pool + energy rank
selects a better (square-planar) conformer. `--ensemble-size 5` fixed BDNN.

`ensemble_size` is plumbed: `verify_roundtrip.py --ensemble-size` → `OIN3DGenerator(...)` /
`MetalloGenAdapter(ensemble_size=...)` → `generate_3d_structures(..., pool_size / num_conformers)`.

## Task

1. **Confirm robustness** (not a one-off): run BDNN several times at the chosen size.
   ```
   for i in 1 2 3; do
     uv run python tests/integration/verify_roundtrip.py --only BDNN \
       --optimizer mace-omol-0-extra-large-1024 --ensemble-size 5 --output-dir /tmp/bdnn_$i \
       2>/dev/null | grep -E "PASS|FAIL|coord sphere\)"
   done
   ```
   Expect PASS (`[Pd_SPL]`, RMSD < 1.0) each time. If it flakes, bump to 8/10 or go to step 3.
2. **Decide how to land it.** Two options — pick per results:
   - **(a) Preferred workflow (user's stated convention):** keep the default at 1 and run
     roundtrip verification with `--ensemble-size 1` first, then `--ensemble-size 5` only for
     complexes that fail in a way ensemble can fix. If landing as docs, update
     `tests/run_verification.sh` help / README to state this convention. Lowest-risk.
   - **(b) Make `ensemble_size>1` the generator default** (e.g. default 3–5 in
     `MetalloGenAdapter`/`generate_3d_structures`) if the added runtime is acceptable. Note
     this multiplies MACE cost per complex — measure before committing.
3. **If the flag is unreliable** (BDNN still flakes at high ensemble): the real fix is
   **geometry-code-aware selection** — among the FF/MACE pool, prefer the conformer whose
   re-perceived geo-code matches the target (`SPL`), not merely the lowest energy (a
   distorted `TPY` can be energetically competitive for floppy N-donors). Implement in
   `generate_3d_structures` selection: after building the pool, encode each survivor's
   coordination geometry (reuse `OINDiscreteAligner`/`_find_best_geometry_match`) and tie-break
   energy ranking by geo-code match to the requested geometry. This is the Opus-escalation path.

## Verification / done criteria

- BDNN round-trip PASSES reliably (≥3/3 runs) at the chosen ensemble size.
- No regression: a full `tests/run_verification.sh --optimizer mace-omol-0-extra-large-1024
  --ensemble-size 1` (then re-run failures at `5`) — roundtrip ≥ 19/25, and if BDNN is
  landed as default-ensemble, ≥ 20/25.
- Unit suites green: `uv run python -m unittest discover tests` (55) and `discover tests/unit`
  (127, skip 3). Encoder `verify_xyz_to_oin.py` still 25/25.
- Log the outcome to `spec/worklog/NOTES.md` and commit `--no-verify` (no push).

## Guardrails

- Round-trip run convention (user): `--ensemble-size 1` first, then `5` only to rescue an
  ensemble-fixable failure. MACE runs are slow (~minutes/complex); use `--only <name>` for
  targeted runs, omit `--optimizer` for a fast FF string-only check.
- Do not touch the string-fidelity code (stereo carry, `_smilesAtomOutputOrder`) — those are
  done and committed. This task is geometry-selection only.
