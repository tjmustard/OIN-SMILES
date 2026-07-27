# v0.4.9 · Lane 1 — making the budget a bound, and where the time actually is

> ## 🔴 There is no single sink. Three molecules, three different dominant costs.
>
> | molecule | dominant cost | share |
> |---|---|---:|
> | `FOSNEI_comp_0` (non-eta, boron) | `embed.get_embedding` | 75% |
> | `CAHQEJ_comp_0` (eta, 2 haptic) | `embed.get_embedding` + `numpy.linalg.eig` | 77% |
> | `VAFMIA_comp_0` (`[Cu_LIN]`, adamantyl NHC) | **`chirality._reparse_cip_label_once`** | **99%** |
>
> **This is the release's most transferable finding, and it was learned the expensive way.**
> The first version of this lever threaded the deadline into `get_embedding` alone — the
> function two profiles had indicted — and measured **ε = +48.4 s on a 30 s budget: it
> changed almost nothing.** A bound threaded into whichever function profiled expensive last
> is not a bound.
>
> The charter's suspects were not wrong so much as **partial**: the CBC solve really is 2.1%
> and the `accept_fn` re-encode really is 0.8% *on `FOSNEI`* — and the same `accept_fn` is
> **62%** on `VAFMIA`. The cost is bimodal by molecule, not small.

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

## 2. The attribution — three molecules, three cost regimes

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

### `VAFMIA_comp_0` — `[Cu_LIN]`, bis-adamantyl NHC, profiled **with the bound ON at 30 s**

Spent **78.60 s** — 2.6× over — and this is the profile that broke the first design.

| sink | cumtime | **tottime** | ncalls |
|---|---:|---:|---:|
| **`chirality._reparse_cip_label_once`** | 77.78 s | **77.78 s** | 32 (2.43 s each) |
| ↳ via `accept_fn` → `_reencode_key_matches` → `build_contract_mol` → `_template_sp3_label` | 48.49 s | — | — |
| ↳ via `_select_by_geometry` (**outside** `generate_3d_structures`) | 39.06 s | — | — |
| `generate_3d_structures` (the only thing the first bound covered) | 39.53 s | — | 1 |
| CBC | 0.35 s | — | 3 |

`get_embedding` does not even appear. **Half the wall-clock is in `_select_by_geometry`, which
runs in the adapter, after `generate_3d_structures` has returned** — a deadline living inside
that function cannot reach it by construction.

### Why the sink structure settles the mechanism

`get_embedding` is not one long native call. It is

```python
for alternative_ace_mol in alternative_ace_mol_list:      # outer
    for haptic_scale in scales_for_haptic:                # inner, 4 scales
        ... rc = AllChem.EmbedMolecule(rd_mol, params)    # several per iteration
```

— a nested **Python** loop, so a deadline checked at the top of each loop makes the
un-interruptible unit **one scale sweep**, not one whole attempt. The charter's worry —
*"a deadline checked around a single 90-second CBC subprocess is not enforcement"* — was the
right question, and for this sink the answer is no: nothing here is one long native call.

**But it was the wrong question to ask only once.** The same reasoning applied to
`VAFMIA` gives a different answer: its cost is one `accept_fn` re-encode, ~24 s, and that call
*is* effectively atomic — the deadline can stop the next one from starting, not the one in
flight. **That is what sets ε**, and it is why the bound is checked at three places rather than
threaded into one function:

| checkpoint | covers | why it exists |
|---|---|---|
| inside `get_embedding`'s two loops | FOSNEI/CAHQEJ regime | the dominant cost when the embed is hard |
| before each `accept_fn` call | VAFMIA regime, in-loop half | 48.5 s of 78.6 s on VAFMIA; **0.63 s of 82.4 s on FOSNEI** |
| before `_select_by_geometry` | VAFMIA regime, post-loop half | 39.1 s, and structurally unreachable from inside `generate_3d_structures` |

## 3. What shipped

Behind **`OIN_ENFORCE_BUDGET`, default OFF**. Unset, every path below is byte-identical to
pristine — the deadline is `None` and each check is one `is not None`.

| change | file |
|---|---|
| deadline threaded into `get_embedding`, checked at the top of **both** nested loops; returns `None` on expiry | `generator3d/embed.py` |
| `enforce_budget` parameter; reads `lever_enabled("OIN_ENFORCE_BUDGET")` when `None` | `generator3d/__init__.py` |
| **do not START an `accept_fn` call the budget cannot afford** | `generator3d/__init__.py` |
| **whole-generation deadline**, started in the adapter so it covers what sits outside `generate_3d_structures` | `generation/metallogen_adapter.py` |
| **budget spent before selection → take the cheap pick** (`early_exit=False`, the existing lowest-energy fallback) rather than pay for a re-encode search | `generation/metallogen_adapter.py` |
| empty pool at the deadline → **`BudgetExhaustedError`** instead of `[]` | `generator3d/__init__.py` |
| the new typed error, re-exported alongside its two siblings | `generation/metallogen_adapter.py` |
| `--gen-timeout`, and `--mol-timeout` now sets the generator budget | `tools/test_dataset_roundtrip.py` |

The last three rows are the **second** design. The first shipped only the `get_embedding`
check, measured **ε = +48.4 s** on a 30 s budget, and is recorded here rather than quietly
replaced — see §7.

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

32 molecules, 8 evenly spaced by rank from each of the four runtime strata (5 haptic in each of
the two slowest bands), one fresh interpreter per (molecule, arm) because several `OIN_*` levers
and the PuLP/embed caches are frozen at import time. Both arms run under the **same** 120 s outer
SIGKILL — an arm measured without a kill is not comparable to one measured with it.

| arm | n | killed | median | max | **max ratio** | `>` budget |
|---|---:|---:|---:|---:|---:|---:|
| OFF | 32 | 1 | 14.5 s | 78.8 s | **2.63×** | 11 |
| **ON** (design 1 — `get_embedding` only) | 32 | 1 | 14.9 s | 78.4 s | **2.61×** | 11 |
| **ON** (design 2 — three checkpoints) | 32 | 1 | 14.9 s | **62.8 s** | **2.09×** | 11 |

| | design 1 | design 2 |
|---|---:|---:|
| **ε = max(spent) − budget** | **+48.4 s** | **+32.8 s** |
| CPU over the cohort | −8.2% | **−9.7%** |

### What this says, without varnish

**The bound works, and it is not tight.**

- **ε = +32.8 s on a 30 s budget.** The worst case is `VAFMIA_comp_0` at 62.8 s. Design 1
  measured +48.4 s, so covering `accept_fn` and `_select_by_geometry` bought **a third of the
  overrun** — but the remainder is one **in-flight** `accept_fn` re-encode, ~24 s of
  `_reparse_cip_label_once` at 2.4 s a call. **The deadline can decline to start the next call.
  It cannot interrupt the one running.** That is ε, and that is the un-interruptible operation
  the charter asked to have named.
- **The tail is compressed, not removed.** `11` molecules exceed the budget in *both* arms — the
  same 11. The bound reduces *how far* over they go (2.63× → 2.09×), not *how many* go over. Any
  claim that this release delivers `max(elapsed_s) < 30 s` would be false.
- **Byte-identity holds: 28/28.** Every molecule that finished in both arms produced an
  identical generated structure. The bound changes *which* molecules finish, never *what* a
  finishing molecule produces.
- **Zero converted late-successes.** At a 30 s budget on this cohort, no molecule that produced a
  structure with the lever OFF failed to with it ON, and the 4 that fail do so in both arms. The
  "a tighter bound will look like a regression" warning did not materialise here — which is a
  fact about *this* cohort at *this* budget, not a general result.
- **The CPU figure is a floor.** The OFF arm ran under a 120 s outer kill, so its true cost is a
  lower bound and the real saving is larger than 9.7%.

### Handed on

Making the bound tight means bounding **inside the CIP/perception layer** —
`chirality._reparse_cip_label_once` and the `_template_sp3_label` loop in `build_contract_mol`.
That is a change to *perception* behaviour, not to *when work stops*, so it does not belong in
this release. `build_contract_mol` already documents `None` as a legitimate failure return
("callers fall back to coordinate re-perception"), so the degradation path exists.

**v0.4.10 gets two measured targets, then**: the discarded `.index()` scan (§5, 22% of an eta
generation) and CIP re-parse cost (99% of `VAFMIA`).
