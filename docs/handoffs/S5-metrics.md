# S5 — RMSD mapping sentinels (996/999) + results hygiene

Branch: `feature/roundtrip-metrics` · Read `docs/handoffs/README.md` first.

## Mission

41 registry rows "fail" with RMSD 996.0/999.0 — but those are SENTINELS from
`tests/integration/rmsd_utils.py`, not geometry: 996.0 = coordination-sphere
per-element atom-count mismatch between the two mols (mapping can't even
start), 999.0 = hard failure inside the metric. A structure can be chemically
correct and still hit 996 (known case: `ABETIK` — correct OIN, coord-sphere C
count 12 vs 13 from an eta-perception difference). Success = the metric
distinguishes "genuinely bad geometry" from "mapping infeasible", sentinel rows
get an honest class, and borderline cases get adjudicated.

## Evidence pack

Fresh sentinel rows: `DAPZIF_comp_0` [Pd_SPY] 996, `CAWYOR_comp_0` [Mn_SPY]
996, `RUBNUZ_comp_0` [Ni_SPL] 996, `ROJXIY_comp_0` [Mn_TPY] 996,
`YICXEO_comp_0` [Y_SPY] 999, `NEZWEU_comp_0` [Sc_OCT] 999.
Stale-known: phospholyl/η⁵ metallocenes `ADALUM`, `ADAMAT`, `AFIHEC`, and
`ABETIK` (Zr Cp/allyl/borate). Y-metal complexes fail 96% overall (27/28) —
eta-heavy chemistry; expect several to be this class or S2's.
Borderline: `SIMNUZ_comp_0` [Ru_TET] RMSD 1.08 — just over the 1.0 threshold;
adjudicate (eta centroid handling? genuinely floppy ligand? threshold?).
Triage bucket also assigned to you: `geometry_or_fragment_change` (15 rows —
generated structure re-encoded with a different metal geometry code or
fragment count; likely a ligand detached during generation → decide whether
that's a real generation defect class for a future wave or a perception issue).

## Also in scope: hard per-molecule watchdog (the `--mol-timeout` is not enough)

`--mol-timeout` uses `signal.alarm` (SIGALRM) + a Python handler. That CANNOT
interrupt a hang inside native C++ (MetalloGen/RDKit conformer search): the
signal is queued until control returns to the Python interpreter, which never
happens. Observed live: `UGUHAH_comp_0` (97-atom, `photo/`) wedged one Phase-0
shard for 35+ min despite a 420 s cap; only an OS `kill` cleared it. Add a real
watchdog — run each molecule in a subprocess and SIGKILL it on wall-clock
timeout (or a `multiprocessing` worker with `.terminate()`), so a single
pathological molecule can't stall a whole run. UGUHAH itself is also a generator
robustness bug worth filing (which native call spins — needs `py-spy`/gdb with
ptrace permission).

## Where the logic lives

`tests/integration/rmsd_utils.py::calculate_tmc_rmsd` — coordination-sphere
RMSD with eta rings reduced to centroids; returns 999.0 on exception (line
~59), 996.0 when the per-element counts of the two coordination spheres differ
(line ~109). The harness treats any value ≥ 1.0 as failure, so sentinels read
as astronomically bad geometry and the report line says "High RMSD: 996.0".

## Design guidance

1. **Honest taxonomy first**: make the harness/report distinguish
   `rmsd_mapping_failed` (sentinel) from `high_rmsd` (real number ≥ 1.0). You
   own `tools/*` — `test_dataset_roundtrip.py` currently formats the error as
   "High RMSD at TIER: 996.0000"; give sentinels their own error string so
   `tools/classify_failures.py` (also yours) can bin them directly.
2. **Then reduce mapping failures**: the per-element count mismatch usually
   comes from eta-ring membership perception differing between the input mol
   (XYZ-perceived) and the generated bonded mol (e.g. η⁵ counted as 5 C on one
   side, 6 on the other). Options: compare on centroid-reduced spheres with
   element-agnostic eta slots; fall back to a metal-anchored assignment
   (Hungarian on distances) when counts differ by ≤2; only sentinel when truly
   unmatchable. Keep it MEAN RMSD (project convention — never max-per-atom).
3. Re-adjudicate the sentinel rows with the improved metric: how many are
   actually fine geometrically (like ABETIK)?

## Verify-first steps

1. Repro `DAPZIF_comp_0` with `--only`; capture both coordination spheres
   (input-perceived vs generated-bonded) and print the per-element counts that
   diverge.
2. Check `YICXEO`/`NEZWEU` 999s: what exception is being swallowed?
3. Only then decide how far to take the mapping fix vs. the honest-label fix.

## Files

- **Own:** `tests/integration/rmsd_utils.py`, `tools/test_dataset_roundtrip.py`,
  `tools/classify_failures.py`, `tools/rebuild_summary.py` (+ new test file
  `tests/unit/test_rmsd_mapping.py`).
- **Read-only:** everything under `src/oinsmiles/` (all other sessions).
- Note: Phase-0 workers and a continuous runner execute
  `tools/test_dataset_roundtrip.py` from the MAIN checkout — your edits live in
  your worktree until merge, so no interference; but rebase before opening the
  PR since Phase 0 already touched these files (commit `d950f2a`).

## Acceptance

- Sentinel rows re-classify as `rmsd_mapping_failed` (or pass with a real
  RMSD after the mapping fix); zero rows report "High RMSD: 996/999".
- `ABETIK_comp_0` adjudicated: either passes with a real RMSD or is documented
  as an eta-perception count issue with the exact atom lists.
- `SIMNUZ_comp_0` adjudicated with evidence.
- Guard tests for: sentinel taxonomy, count-mismatch fallback mapping, eta
  centroid reduction. Full unit suite green.
