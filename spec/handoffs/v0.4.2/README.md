# v0.4.2 — Round-Trip Accuracy Fix Wave · Parallel Session Protocol

Seven fix phases (S1, S3, S5, S6a, S6b, S7, docs) plus a baseline (P0) and a validation
capstone, each in its own git worktree, each owning a **function-granular** slice of the code.
This wave attacks the round-trip **accuracy** defect classes that the full-dataset (>17,000-molecule)
tmCAT/tmPHOTO sweep re-surfaced — the same classes the v0.3.6 (S1–S6) and v0.3.7 (R1–R5) waves
named, now at ~6× scale with their deferred residuals.

It is a **separate stream** from the v0.4.0-perf speed wave (`spec/handoffs/v0.4.0-perf/`). That
wave attacked cost centers; this one attacks encoder/generator accuracy.

**Read your own phase doc, then this protocol. Do not touch files (or functions) another phase
owns** — the ownership matrix below is enforced at function granularity, exactly as v0.3.6's S1/S2
and v0.4.0's P1/P3/P5/P11 were.

Each phase doc is **self-launching**: a `▶ START HERE` bootstrap at its top creates the worktree,
lists what to read, and gives a repro command. Open it in a fresh Claude Code session in the main
checkout and follow it.

**These handoff docs live under `spec/handoffs/` which is GITIGNORED in the main checkout, and are
force-committed onto `release/v0.4.2`** so they also appear inside any worktree branched from the
release branch. Read them by absolute path from the main checkout
`/home/tjmustard/Documents/GitHub/OIN-SMILES/spec/handoffs/v0.4.2/`.

## Staging model — read this before branching

**This wave stages on a dedicated integration branch, `release/v0.4.2` — not on `main`.** It is
based on **`c7edeeb6` (= tag `v0.4.1`, current `main`)**. Every phase **branches from and
squash-merges back into `release/v0.4.2`**, and `main` sees v0.4.2 exactly once — as a single
reviewed squash after the capstone signs off (project convention; PR #2/#3, v0.4.0 precedent).
`release/v0.4.2` and `main` stay **unpushed** unless the user says otherwise (standing instruction).

The integration/staging worktree is **`../OIN-SMILES-v0.4.2`** (on `release/v0.4.2`) — it is the
**merge target** every phase squash-merges into, and the capstone's validation ground. It already
exists.

## The baseline commit — one commit, no exceptions

**`c7edeeb6` is THE baseline for the whole wave.** The v0.4.0 wave's hardest-won lesson: a
mixed-provenance floor is not a floor. Do **not** quote any headline pass-percentage from
`results-v0.4.0/` — that accumulator is mixed-provenance (`5538b722-dirty` **and** `c7edeeb6`
stamps) and `--quick`. The only valid floor is `spec/handoffs/v0.4.2/BASELINE.md` (committed by
P0): a **set of molecule IDs that pass on `c7edeeb6`**, plus per-class goldens. Gate per-molecule
against that set, never on a percentage.

## The live accumulator keeps running (do not fight it)

A `tools/test_dataset_roundtrip.py --quick --continue --random` runner is live, writing to
`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0/` — it keeps **growing** dataset coverage and is the live
backlog source (regenerate `CASE_REGISTRY.md` / `V0.4.1_ACCURACY_BACKLOG.md` any time with
`tools/classify_failures.py` + `tools/group_v041_backlog.py` **on a copy**). It is **not** the
floor.

### ⚠ Only ONE dataset sweep may run at a time — across ALL worktrees

Concurrent sweeps fabricate `no_conformers` failures and corrupt every timing number (v0.3.6
invented one for `BUYNIU`). The live `--quick` accumulator **counts as a sweep**. Therefore:

- The two **large** sweeps — **P0** (clean floor) and the **capstone** (integrated A/B) — must
  **pause the accumulator** while they run, then resume it. (Find it with
  `pgrep -af test_dataset_roundtrip`; pause = stop the `--continue` loop; resume = restart it.)
- **Per-phase** regression checks use **targeted `--only <goldens>` repros into a private `/tmp`
  dir + the unit suite + byte-identity** — NOT full sweeps. This is the v0.3.6/v0.3.7 acceptance
  model ("re-run the named currently-passing molecules"). A single flaky `--only` row is
  re-verified **alone** before it is believed.
- **`classify_failures.py --output-dir` WRITES to the dir you give it.** Never aim it at
  `results-v0.4.0/`; run it on a copy.

## Worktree setup

```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES
git worktree add ../OIN-SMILES-<slug> -b feature/roundtrip-<slug> release/v0.4.2
cd ../OIN-SMILES-<slug>
uv sync
```

`uv sync` in a fresh worktree installs **rdkit 2026.3.3**; the long-running baseline venv
(`~/.venv`, used by `uv run --no-sync` in the main checkout) is **2025.09.3**. Both are blessed;
the golden-env guard test accepts the set. **The two disagree on `/`\`` double-bond direction** for
some `smiles_1` (e.g. DIXJOL/NOLPOV) while being canonical-key-equal — so run any `smiles_1` diff
under **one** pinned rdkit or you will chase phantom drift (S3 bisected six commits before spotting
this). **One clone, many worktrees — do NOT make a second clone.**

## Reproducing a case

The dataset lives only in the main checkout (gitignored). From your worktree:

```bash
uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir /home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset \
    --output-dir /tmp/rt-<slug> --quick --only <MOLECULE>[,<MOLECULE>...]
```

The report JSON `<output-dir>/individual_reports/<MOL>.json` carries `smiles_1` (input encode),
`smiles_2` (re-encode of generated 3D), the error, and provenance (`commit_id`, `rdkit_version`).
**Never point `--output-dir` at `results-v0.4.0/`.** `--only` matches on basename, and some
molecules exist in both `photo/` and `cat/` — reports are keyed by basename, so the two overwrite
each other.

**`--quick` measures a different generator** (`ff_params_fast = {"uff_pool_size": 2,
"max_attempts": 10}`, 30 s hard-kill). It is fine for a fast string/atom-count repro of a
*deterministic* defect, but for anything RMSD- or conformer-yield-sensitive (S5 winding, S7
metrics, no_conformers) re-run **non-`--quick`** with `--mol-timeout ≥ 1700s` before believing a
result. A genuine `no_conformers` row (`ZIHGEE`) burns ~1696 s; a `--mol-timeout` under ~600 s
SIGKILLs healthy-but-slow molecules and you blame your diff.

## The two cost-free correctness proofs

- **Byte-identity to pristine** (for pure-refactor hunks): generation is seeded (`seed=42` default,
  deterministic since v0.4.0), so a refactor that changes no chemistry produces a **bit-for-bit
  identical** re-encode/XYZ to pristine at the same seed. Stronger and cheaper than any sweep. Use
  `git show HEAD~1:<path> > <path>` … `git checkout -- <path>` to A/B — **never** `git stash push`
  an already-committed (clean) path (it saves nothing, and `stash pop` then applies the long-lived
  `!!GitHub_Desktop<…>` stash into your tree — this has bitten three sessions).
- **Simulate the contract mol instead of generating one**: `Chem.Kekulize(tmc_mol)` (flags
  retained) is exactly what `build_contract_mol` hands `get_oin_string`. That turns a multi-minute
  round trip into a ~1 s check against committed fixtures — how S3 guard-tested its core defect
  without touching the dataset.

## Phase / file-ownership matrix (function-granular)

| Phase | slug | Attacks | Owns (edit freely) |
|---|---|---|---|
| P0 | `baseline` | clean floor + goldens | `spec/handoffs/v0.4.2/BASELINE.md` — **no src** |
| S1 | `donor-h` | donor_H_atom_count (fixable subset), H_on_terminal_oxo_imido | `oin/inline.py` (`:265-317`); `metallogen_adapter.py::convert_parsed_to_msmiles` H-strip `:164-222` (incl. the `:114-116` NON raise); `generator3d/ligand.py::get_ligand_from_smiles` — **the `AddHs` line `:63`** |
| S3 | `aromatic` | encode_crash_other, kekulize_encode_crash, macrocycle_perception, garbled_aromatic | `utils/xyz2mol.py`, `utils/aromaticity.py` — **fully disjoint** |
| S5 | `geometry` | geometry_or_fragment_change, geometry_NON, winding_flip | `utils/oin_aligner.py` geometry matcher/classifier/templates; `metallogen_adapter.py::_select_by_geometry :1121-1213` (**incl. the winding block `:1173-1197`**) + `generate` pool `:1264-1275` + `_eta_winding_multiset :1043` + `OIN_TO_METALLOGEN_GEO` dict `:73-89` (**dict only — NOT the `:114-116` raise, S1's function**) |
| S6a | `ez-stereo` | EZ_bond_stereo | `generator3d/embed.py::_apply_double_bond_stereo :417-469`; `generator3d/ligand.py::get_ligand_from_smiles` — **the `near_donor` block `:55-62` only**; `core/translator.py::_clear_chelate_locked_bond_stereo` (verify) |
| S6b | `atom-stereo` | atom_stereo, `[S@SP3]` string_mismatch subset | `core/chirality.py`; `generator3d/embed.py::_apply_atom_chirality :480 / _permutation_is_odd :470`; `metallogen_adapter.py` sp3 CIP stamp `:559,:776-778`; (audit-only) `ligand.py` chiral-centre capture `:73-83` |
| S7 | `metrics` | high_rmsd + no_conformers/timeout/gen_exception **triage** | `generator3d/bond_lengths.py::ENABLED_METALS`, `generator3d/clean_geometry.py`, `generator3d/embed.py::get_embedding :598` **and the stereo call-sites `:669,:764`**; honesty review in `tools/test_dataset_roundtrip.py` |
| docs | `docs` | carborane_unsupported + notation/artifact residuals | `docs/KNOWN_LIMITATIONS.md`, `spec/handoffs/v0.4.2/wontfix-carboranes.md` |
| capstone | `validation` | integrated A/B; one squash → main | `spec/handoffs/v0.4.2/VALIDATION.md` — **no src** |

Everyone may **ADD** new test files under `tests/unit/` (name them distinctly).

### Shared-function serializations (the only real hazards)

| Shared function | Co-editors | Rule |
|---|---|---|
| `generator3d/ligand.py::get_ligand_from_smiles` | S1 (`AddHs :63`), S6a (`near_donor :55-62`), S6b (`chiral_centers :73-83`, audit) | **S1 lands first**; S6a rebases and edits only `:55-62`; if S6b needs `:73-83`, it chains after S1 too. The captures are index-coupled to the post-`AddHs` frame. |
| `metallogen_adapter.py::_select_by_geometry` (+ `generate` pool) | winding ∩ geometry | **Merged into S5** — single owner. The winding pick at `:1173-1197` consumes S5's `scored` list at `:1180`; they are data-coupled, not disjoint. |
| `generator3d/embed.py` stereo apply-fns | S6a (`_apply_double_bond_stereo`), S6b (`_apply_atom_chirality`), S7 (`get_embedding` + call-sites) | Function-disjoint bodies. **S7 owns the call-sites `:669,:764`**; S6a/S6b edit only the `_apply_*` bodies. |

## Ordering

```
P0 baseline (pin c7edeeb6; goldens; no src) ── land first; everything gates on it
      │
      ├─ S3 aromatic      (xyz2mol/aromaticity — fully disjoint)      ┐ run concurrently,
      ├─ docs             (KNOWN_LIMITATIONS — disjoint)              │ land in any order,
      ├─ S5 geometry      (oin_aligner + _select_by_geometry)         │ rebase on each squash
      ├─ S6b atom-stereo  (chirality + embed _apply_atom_chirality)   │
      ├─ S7 metrics       (bond_lengths/clean_geometry + embed        │
      │                    call-sites — coordinate w/ S6a/S6b)        ┘
      │
      └─ ligand.py serial chain:  S1 donor-h ──► S6a ez-stereo
                                  (owns AddHs)   (rebases; narrows near_donor
                                                  + fixes _apply_double_bond_stereo)
      │
      ▼
   capstone (integrated per-molecule A/B on release/v0.4.2) ──► one squash → main
```

- **After any phase squash-merges into `release/v0.4.2`, live worktrees `git rebase
  release/v0.4.2`** before landing theirs.
- **Tag `archive/<slug>` before any `git branch -D`** — a squash merge means `-D` throws away the
  granular history (S3 and S4 both lost theirs in prior waves). Verify **content**
  (`git diff --quiet release/v0.4.2..feature/roundtrip-<slug> -- <file>`), not ancestry.

## Acceptance criteria (every fix phase)

1. **Sub-triage first, then scope honestly.** Reproduce your class's defects on `c7edeeb6` from
   your goldens; verify which are real accuracy bugs vs notation/harness limits (S1's +3 rows are
   partly a notation limit; S7's high_rmsd/timeout are largely artifacts). **Promise the fixable
   subset, route the rest to `docs`.** Treat your handoff hypothesis as a hypothesis — five prior
   sessions shipped a wrong root cause with the refuting evidence already in the repo.
2. **Your goldens round-trip** (canonical key + atom-count/RMSD where your class applies) via
   `--only` into a private dir; **currently-passing spot-check** (your named passers still pass).
   For pure-refactor hunks, prove **byte-identity to pristine** at a fixed seed.
3. **New unit guard tests** under `tests/unit/`, each failing against pre-fix code.
4. **Full suite green**: `uv run python -m unittest discover tests/unit`. **Measure the baseline
   count in your OWN worktree/venv first** — every count quoted in an older handoff has been stale.
   `ruff check` clean (note: `ruff --fix` on D-rules can truncate docstrings — re-read after
   autofix; the pre-commit hook enforces `ruff format`).
5. **Before/after evidence** (named molecule flips, per-class, RMSD) in the squash-commit body.

## Landing a phase

From the finished phase session (see `SESSION_PROMPTS.md` for the exact prompt):

1. Commit on `feature/roundtrip-<slug>`; `git rebase release/v0.4.2` (absorb prior squashes);
   re-run your gate + `unittest discover tests/unit`.
2. `cd ../OIN-SMILES-v0.4.2 && git merge --squash feature/roundtrip-<slug> && git commit`
   (subject `v0.4.2-<phase> <slug>: <one-line>`; body = named flips + regression result).
3. Confirm content landed: `git diff --quiet release/v0.4.2..feature/roundtrip-<slug> -- <owned files>`.
4. `git tag archive/roundtrip-<slug> feature/roundtrip-<slug>`; `git worktree remove ../OIN-SMILES-<slug>`; `git branch -D feature/roundtrip-<slug>`.
5. Announce so live phases `git rebase release/v0.4.2`.

## Environment gotchas (carried from v0.3.6/v0.4.0 — all still true)

- The metal is **always `fragments[0]`** in an OIN (load-bearing invariant).
- Use **mean** coordination-sphere RMSD, never max-per-atom, for any metric threshold. RMSD
  sentinels ≥900 mean "metric could not run", not "bad geometry".
- No `xtb` binary → `optimizer="g-xtb"` (the default) warns and returns FF geometry unchanged. This
  wave is **FF-only by decision**; do not chase optimizer-level RMSD noise. Construct generators
  with `optimizer=None` when benchmarking.
- MetalloGen's geometry key `8_squre_antiprismatic` is **misspelled upstream** — do not "correct"
  it (`metallogen_adapter.py:86-91` byte-matches it deliberately; it broke once).
- `MolToSmiles` ignores stored `SetNumExplicitHs` on an inherited atom (put H in isotope);
  `SanitizeMol` rejects 4-coordinate neutral boron; `FastFindRings` is **not** SSSR.
- `status: pending_g-xtb` is **not** a failure; conformer selection is stochastic across *unseeded*
  runs but the default `seed=42` is deterministic — a fixed-seed flip is a real chemistry change, an
  unseeded one is not.
