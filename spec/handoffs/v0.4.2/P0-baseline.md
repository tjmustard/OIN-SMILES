# ▶ START HERE — P0 baseline (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** P0 is the
**prerequisite for every other phase** — it produces the clean floor and the per-class goldens all
phases A/B against. It touches **no source**, only `BASELINE.md`.

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-baseline -b feature/roundtrip-baseline release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-baseline && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `/home/tjmustard/Documents/GitHub/OIN-SMILES/spec/handoffs/v0.4.2/README.md`
- the live backlog — `.../tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0/V0.4.1_ACCURACY_BACKLOG.md`, `CASE_REGISTRY.md`
- the harness — `tools/test_dataset_roundtrip.py` (gates: string → atom-count `:177-184` → RMSD `:173`; `--quick` budget `:616-617`), `tools/classify_failures.py` (`SESSION_OF`, `--output-dir` WRITES)

### 3 · Mission

Produce `spec/handoffs/v0.4.2/BASELINE.md`: a **single-commit** (`c7edeeb6`) floor and per-class
goldens. **No headline percentage** — the deliverable is a *set of molecule IDs* + goldens.

1. **Pause the live `--quick` accumulator** (see `SESSION_PROMPTS.md` → "Pausing / resuming"). Only
   one sweep at a time.
2. **Re-sweep only the fixable failing classes on `c7edeeb6`**, serially, **non-`--quick`**, with
   `--mol-timeout 1800`:
   - classes: `donor_H_atom_count`, `H_on_terminal_oxo_imido`, `geometry_NON`,
     `geometry_or_fragment_change`, `winding_flip`, `EZ_bond_stereo`, `atom_stereo`,
     `encode_crash_other`, `kekulize_encode_crash`, `macrocycle_perception`, `garbled_aromatic`,
     and the `[S@SP3]` subset of `string_mismatch_other`.
   - Pull the current member list per class from a **fresh** `classify_failures.py` run on a **copy**
     of `results-v0.4.0/`, then feed the union as `--only A,B,C,...` into a **private**
     `--output-dir /tmp/v042-baseline` (never the shared dir).
   - `--quick` measures a different generator (`uff_pool_size=2, max_attempts=10`, 30 s) — do not use
     it for the floor.
3. **Do NOT re-sweep the pure-artifact classes for the floor** — `timeout` (806), `high_rmsd` (94),
   `carborane_unsupported` (79), `no_conformers` (225). Highest cost, lowest diagnostic value; they
   belong to S7's own full-fidelity triage and to `docs`. (Re-sweeping 806 timeouts alone is ~13 h;
   genuine `no_conformers` rows push it to days.)
4. **Provenance-filter the passing set**: from `results-v0.4.0/`, the "must-not-regress" set is the
   molecules already stamped `commit_id == c7edeeb6` **and** `status == success`. Flag
   `5538b722-dirty` passers as **untrusted** (a molecule that "passed" on the dirty tree may not on
   `c7edeeb6`) and list them separately.
5. **Freeze per-class goldens**: for each fixable class, pick 4–8 representative molecules (span
   metals/CN/size; include the backlog's named representative) that reproduce the defect on
   `c7edeeb6`. These become each phase's A/B set.
6. **Resume the accumulator** exactly as it was.

### 4 · `BASELINE.md` must contain
- The baseline commit (`c7edeeb6`), rdkit version used, date, and the exact sweep command.
- Per fixable class: current count, the frozen goldens (with the one-line defect signature), and the
  repro'd-on-`c7edeeb6` confirmation.
- The **must-not-regress passing set** (IDs stamped `c7edeeb6`), and the **untrusted-provenance**
  list, explicitly separated.
- An explicit "**no headline % here**" note and a pointer to why (README → "The baseline commit").
- The artifact-class counts (timeout/high_rmsd/carborane/no_conformers) as *context only*, marked
  "not part of the floor; S7/docs own them."

### 5 · Gate & landing
- Gate: `BASELINE.md` exists, single-commit, with goldens that reproduce on `c7edeeb6` and a
  provenance-filtered passing set. No src touched.
- Optionally (user's call) kick off a **background full 14,771+ sweep on `c7edeeb6`** as capstone
  insurance — but do **not** block the wave on it.
- Land: squash-merge `BASELINE.md` into `release/v0.4.2` (see `SESSION_PROMPTS.md`). It lands first;
  everything else gates on it.
