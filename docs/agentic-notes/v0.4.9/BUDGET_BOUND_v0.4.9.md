# v0.4.9 · Lane 1 — making the budget a bound, and where the time actually is

> **The charter named two suspects. Profiling indicts neither.**
>
> The unbounded PuLP/CBC solve is **2.1%** of a generation. The 48–57 s `accept_fn`
> re-encode is **0.8%**. The sink is `embed.get_embedding` — **61.5 s of *self* time out of
> 82.4 s** — and because it is a nested **Python** loop over `AllChem.EmbedMolecule` rather
> than one long native call, a deadline checked inside it is real enforcement. No
> `fork`/`RLIMIT_CPU` machinery is needed, and none was added.

## 1. What was wrong, precisely

`OIN3DGenerator(timeout=)` is handed down as `embed_time_budget` and becomes a deadline
that is checked **only at the top of the embed attempt loop**
(`generator3d/__init__.py:350`, enforced `:457`/`:529`). An in-flight attempt always runs to
completion. Two independent probes measured the consequence, both without an outer kill:

| source | requested | observed | worst |
|---|---:|---|---:|
| eta sample | 60 s | 60.7 – 137.9 s | `GOHWOQ` **2.3×** |
| boron set, 33 molecules | 60 s | 60.0 – 172.8 s | **2.9×** |

**What is *not* evidence for this**, despite being the release's chartered justification:
`max(elapsed_s) = 759.9 s`. That field is a sum over up to three separately SIGKILLed
harness attempts, and all 4658 single-attempt rows in the 5k sweep finish within **0.2 s**
of their 300 s cap. Full refutation: `ELAPSED_S_IS_A_SUM_v0.4.9.md`.

## 2. The attribution — two molecules, two populations

`tools/profile_eta.py`, replicating the harness `UFF_1` tier (`optimizer=None`,
`ensemble_size=1`) at a 300 s budget, on a quiet box.

### `FOSNEI_comp_0` — non-eta, boron cage, the corpus's 759.9 s worst case

Generated **successfully in 82.44 s**.

| sink | cumtime | **tottime** | ncalls |
|---|---:|---:|---:|
| `generate_3d_structures` | 82.14 s | — | 1 |
| **`embed.get_embedding`** | 64.83 s | **61.51 s** | 10 |
| `clean_geometry.ff_clean` | 16.69 s | **15.18 s** | 10 |
| CBC (`solve_CBC` / `actualSolve`) | 1.74 s | — | 18 |
| `_reencode_key_matches` (`accept_fn`) | 0.63 s | — | 15 |

### `CAHQEJ_comp_0` — eta, `[Ni_TPL]`, two haptic ligands, 300.2 s in the sweep

Generated **successfully in 72.02 s**.

| sink | cumtime | **tottime** | ncalls |
|---|---:|---:|---:|
| **`embed.get_embedding`** | 55.16 s (77%) | **25.70 s** | 34 |
| `numpy.linalg.eig` via `chem.get_c_eig_list` | 15.99 s | **15.11 s** | 198 |
| `embed._finalize_positions` | 10.35 s | — | 90 |

**Both populations agree on the sink.** CBC and the `accept_fn` re-encode are noise in both.
The "48–57 s per re-encode call" figure that made `accept_fn` a suspect does not generalise —
15 calls here cost 0.63 s in total.

### Why that settles the mechanism

`get_embedding` is not one long native call. It is

```python
for alternative_ace_mol in alternative_ace_mol_list:      # outer
    for haptic_scale in scales_for_haptic:                # inner, 4 scales
        ... rc = AllChem.EmbedMolecule(rd_mol, params)    # several per iteration
```

— a nested **Python** loop. A deadline checked at the top of each loop makes the longest
un-interruptible unit **one scale sweep**, not one `get_embedding` (6.5 s mean on FOSNEI) and
certainly not one whole attempt. The charter's worry — *"a deadline checked around a single
90-second CBC subprocess is not enforcement, it is the current advisory behaviour with more
code"* — was the right question to ask and the answer is no: nothing here is a 90-second
un-interruptible call.

## 3. What shipped

Behind **`OIN_ENFORCE_BUDGET`, default OFF**. Unset, every path below is byte-identical to
pristine — the deadline is `None` and each check is one `is not None`.

| change | file |
|---|---|
| deadline threaded into `get_embedding`, checked at the top of **both** nested loops; returns `None` on expiry | `generator3d/embed.py` |
| `enforce_budget` parameter; reads `lever_enabled("OIN_ENFORCE_BUDGET")` when `None` | `generator3d/__init__.py` |
| empty pool at the deadline → **`BudgetExhaustedError`** instead of `[]` | `generator3d/__init__.py` |
| the new typed error, re-exported alongside its two siblings | `generation/metallogen_adapter.py` |
| `--gen-timeout`, and `--mol-timeout` now sets the generator budget | `tools/test_dataset_roundtrip.py` |

Three design choices worth stating:

- **Expiry returns `None`, it does not raise, inside `get_embedding`.** An exhausted budget
  then looks exactly like an embed that produced no positions, so the attempt loop's
  bookkeeping is untouched and *one* place decides what it means.
- **A non-empty pool is still returned.** A bound should stop work, not discard an answer it
  already has.
- **`StructuralAssemblyError` outranks `BudgetExhaustedError`.** A uniformly structural
  failure would have happened at any budget, so it is the better diagnosis; relabelling it
  "out of time" would send the next release chasing a compute problem that does not exist.

`BudgetExhaustedError` exists so v0.4.10 can tell **its own regressions** from **this
release's intended behaviour**. Without it both arrive as `MetalloGen failed to generate any
conformers`, which `tools/classify_failures.py` buckets as `no_conformers`.

## 4. The plumbing defect that made the A/B impossible

`tools/test_dataset_roundtrip.py:808` hardcoded `timeout_val = 30 if args.quick else 300` and
passed *that* to `OIN3DGenerator(timeout=)`. `--mol-timeout` reached only `_supervise()`'s
SIGKILL. **The two budgets were fully decoupled**: asking for `--mol-timeout 60` still handed
the generator 300 s, so the kill always won and the generator's own budget was never the thing
under test.

Now `--mol-timeout` sets both, `--gen-timeout` overrides, and both are stamped into every
report. Behaviour-identical for every invocation this project has actually run — `--quick`
→ 30, `--mol-timeout 300` → 300, no flag → 300 — and only `--mol-timeout N` for N ≠ 300
changes, which is the fix.

## 5. Handed to v0.4.10: a measured, one-line, 22% win

Not landed here. v0.4.9 changes *when work stops*, not *how long it takes*, and landing an
optimization inside the release that claims to optimize nothing would muddy its own A/B.

`generator3d/embed.py`, top of `get_embedding`'s outer loop:

```python
for alternative_ace_mol in alternative_ace_mol_list:
    alternative_ace_mol_list.index(alternative_ace_mol)   # <-- result DISCARDED
```

`list.index` does a linear scan calling `Molecule.__eq__` → `is_same_molecule` →
`get_c_eig_list` → `numpy.linalg.eig`. On `CAHQEJ_comp_0` that is **3711 `index` calls, 198
eigendecompositions, 15.99 s of a 72.02 s generation — 22%, on a value that is thrown away.**

`is_same_molecule` is pure (it compares eigenvalue lists and returns a bool), so deleting the
line is byte-identical by construction. **Verify that claim before landing it** — this project
has been bitten by "obviously inert" changes before.

## 6. Reproducing

```bash
cd <checkout>
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python   # rdkit pinned ==2025.9.3
export PYTHONPATH=$PWD/src
D=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset

# §2 -- the attribution (quiet box; these are ~80 s each)
$V tools/profile_eta.py $D/cat/FOSNEI_comp_0.xyz --timeout 300 --top 40
$V tools/profile_eta.py $D/cat/CAHQEJ_comp_0.xyz --timeout 300 --top 30

# §7 -- the A/B: asked for B, how many seconds were spent?
$V tools/budget_bound_ab.py --cohort-dir $D/cohort-v049-strata \
    --names-file <32 stratified names> --budget 30 --off-arm-kill 120 \
    --out $D/results-v0.4.9-bound/ab_budget30.jsonl

$V -m unittest tests.unit.test_budget_bound -v      # 10 tests
$V -m unittest discover tests/unit -q               # 930 OK
```

## 7. The A/B — asked for 30 s, how long did it take?

32 molecules, 8 evenly spaced by rank from each of the four runtime strata (5 haptic in each
of the two slowest bands), one fresh interpreter per (molecule, arm) because several `OIN_*`
levers and the PuLP/embed caches are frozen at import time. Both arms run under the **same**
120 s outer SIGKILL, so they are measured under identical outer conditions — an arm measured
without a kill is not comparable to one measured with it.

*(Results table, ε, and the byte-identity check: see `LANE-budget-bound.md`.)*
