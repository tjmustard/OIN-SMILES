# ▶ START HERE — S7 metrics + robustness triage (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** S7 is a
**triage + honesty** phase, not a promise to clear 806 timeouts and 94 high_rmsd — most of those are
**harness/FF-only artifacts, not chemistry**. Deliver the one real win (`ENABLED_METALS`), separate
genuine failures from artifacts, and make the degradation honest.

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-metrics -b feature/roundtrip-metrics release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-metrics && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- prior — `spec/handoffs/v0.4.0-perf/P4-bond-lengths.md` (the per-metal table + `ENABLED_METALS`
  gate), `P6-harness-honesty.md` (**relabel, don't delete**), `P8-clean-geometry-regression.md`,
  `docs/KNOWN_LIMITATIONS.md:154-192` (the 9 genuine `no_conformers`)

### 3 · Verified code paths & what's real vs artifact
- **high_rmsd (94):** HARNESS threshold + FF bias. Gate `tools/test_dataset_roundtrip.py:173`
  (`if rmsd >= 1.0`) runs **after** the string round-trip already matched — so **every high_rmsd case
  is a chemically-correct round-trip** failed on geometric tightness under FF-only (no `xtb`). The
  real, bounded win: metals outside `generator3d/bond_lengths.py:108 ENABLED_METALS =
  {Ni,Pd,Pt,Zn,Cd,Hg,Ag}` fall back to covalent-radius-sum × scale, which **systematically
  overestimates** dative bonds (`bond_lengths.py:5-6`) and inflates coordinate-sphere RMSD.
- **timeout (806):** LARGELY a `--quick` artifact — `test_dataset_roundtrip.py:616-617`
  (`timeout_val=30`, `ff_params_fast={uff_pool_size:2, max_attempts:10}`) + SIGKILL `:253/:302`. A
  full-budget run (300 s / 250 attempts) recovers many. Not a chemistry bug.
- **no_conformers (225):** a **mix** — genuine perception/valence gaps (9 documented,
  `KNOWN_LIMITATIONS.md:154-192`) plus quick-starved embeds. Embed `embed.py:598 get_embedding`
  (retry ladder; `EmbedMolecule` returns **-1** on failure, not an exception — a counter that only
  catches exceptions reports 100% success); surfacing `metallogen_adapter.py:1302-1305`.
- **gen_exception_other (61):** valence/PuLP exceptions + `UncoordinatedFragmentError`
  (`metallogen_adapter.py:92,:154-157` — outer-sphere counterions/solvents, a representation gap).

### 4 · Mission & scope guard (do NOT chase the artifacts)
1. **The real win — extend `ENABLED_METALS`** (P4 precedent): add real per-(metal, ligand) bond
   lengths for the metals that actually appear in high_rmsd goldens but are outside the current set.
   Gate: median coordination-sphere **mean** RMSD improves on those metals, **no per-metal bucket
   regresses** (P4/P8 lesson — a blanket table regressed BINAP; keep it metal-gated).
2. **Honest triage, not a percentage chase.** Re-run a sample of high_rmsd / timeout / no_conformers
   **non-`--quick`, `--mol-timeout 1800`, serially** (pause the accumulator) to separate genuine from
   artifact. Report the split. Do **NOT** loosen the 1.0 Å RMSD gate to "fix" high_rmsd, and do
   **NOT** count `--quick` timeouts as accuracy failures.
3. **Harness honesty** (P6 precedent — relabel, don't delete): make the FF-only degradation and the
   `--quick` budget legible in the harness output/report (e.g. distinguish "FF-floor high_rmsd" and
   "quick-timeout" from real failures), so the backlog stops conflating them.
4. **Genuine no_conformers**: fix the few that are real perception/valence gaps in `get_embedding` if
   tractable; hand `UncoordinatedFragmentError` (outer-sphere) to `docs` as a representation limit.
5. Everything you can't fix → **route to `docs`** with the reason.

### 5 · Owned files (edit only these)
- `src/oinsmiles/generator3d/bond_lengths.py` (`ENABLED_METALS` + table), `generator3d/clean_geometry.py`.
- `src/oinsmiles/generator3d/embed.py` — **`get_embedding :598` and the stereo call-sites `:669,:764`
  only** (you own the call-sites; S6a/S6b own the `_apply_*` bodies — coordinate, don't edit theirs).
- `tools/test_dataset_roundtrip.py` — labeling/honesty only; do **not** change the gate thresholds
  without flagging it as a deliberate harness-correctness change in the commit body.

### 6 · Gate
- `ENABLED_METALS` win: named per-metal median mean-RMSD improvement on the targeted metals, no
  bucket regressing, over N≥10 seeds; the four v0.4.0 goldens byte-identical-or-better.
- Triage report: high_rmsd / timeout / no_conformers each split into genuine vs artifact with counts,
  in the commit body and routed to `docs`.
- New guard tests under `tests/unit/` (bond-length table entries; embed `-1` accounting), failing
  pre-fix. Full unit suite green (own baseline first); `ruff check` clean.

### 7 · Landing
Squash-merge into `release/v0.4.2` (see `SESSION_PROMPTS.md`).
