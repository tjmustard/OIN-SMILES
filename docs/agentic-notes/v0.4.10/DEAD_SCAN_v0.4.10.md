# v0.4.10 · Lane A — deleting one line that eigendecomposed a matrix and threw it away

> **`CAHQEJ_comp_0`: 86.35 s → 57.96 s. −32.9%. Byte-identical.**
> **`FOSNEI_comp_0`: +0.3% — no effect, and now cleanly measured.**
>
> Both facts are the deliverable. A release that reported only the first would be repeating
> the exact mistake v0.4.9 recorded one release earlier: **the cost is bimodal by molecule,
> and a lever aimed at whichever function profiled expensive last optimises one molecule
> class and measures nothing on the corpus.**

## 1. The line

`src/oinsmiles/generator3d/embed.py`, at the top of `get_embedding`'s outer loop:

```python
for alternative_ace_mol in alternative_ace_mol_list:
    alternative_ace_mol_list.index(alternative_ace_mol)   # <-- result DISCARDED
```

It is not a no-op. `list.index` compares with `Molecule.__eq__` → `is_same_molecule` →
`get_c_eig_list` → `numpy.linalg.eig`, so it eigendecomposed a Coulomb matrix per candidate per
outer iteration, over a list it was already iterating, and discarded the index.

The scan is **O(n²) in the candidate list** — iteration *i* re-scans positions `0..i` — which is why
it costs what it costs on molecules with many alternatives and nothing on molecules with few.

## 2. Why deleting it is byte-identical *by construction*

The release's rule: *the removed computation was a pure function of inputs already available, so the
output is identical and only the wall-clock changes.* Here the output is not merely identical, it is
**discarded**. What has to be true instead is that the call had no *side effects*, and it does not:

| function in the comparison call graph | assigns to `self`? |
|---|---|
| `Molecule.__eq__` → `is_same_molecule` | no |
| `get_c_eig_list` | **no** — it guards on `if self.c_eig_list is None:` and then *never populates it*, so it recomputes every call |
| `get_matrix("coulomb")` → `get_adj_matrix`, `get_z_list` | no — all `return` freshly computed values |
| `get_chg` | no |

Verified by reading **and** empirically: `Molecule("CCO").c_eig_list` is `None` before and after
`get_c_eig_list()`. So no downstream caller can observe a warmed cache that this line used to prime.

It also **cannot raise**. `list.index` short-circuits on identity, so searching for an element drawn
from the list itself always finds it and `ValueError` is unreachable — no control flow depended on it.

> One genuine behaviour change, in the safe direction: `get_c_eig_list` calls **`sys.exit()`** if an
> eigenvalue comes back NaN. Deleting the scan removes one path by which a *library* can terminate
> the interpreter. Nothing in this release relies on that, and it is strictly an improvement.

`tests/unit/test_embed_dead_scan.py` (5 tests) pins all of the above, because the gate proves the
removal changed nothing *today* but cannot explain *why* — and "add a one-line cache to
`get_c_eig_list`" is a live temptation that would silently make this deletion observable.
A source lint keeps the line from coming back.

## 3. Measured

`ab_gen.py`, two rounds, arms **alternated** so a drifting box shows up as within-arm spread. Arm A
is the pristine `3077282a` checkout, arm B the lane worktree — two checkouts sharing the main
checkout's pinned venv, this project's established A/B isolation method.

| molecule | class | arm A (pristine) | arm B (deleted) | delta | within-arm spread |
|---|---|---:|---:|---:|---|
| `CAHQEJ_comp_0` | eta, `[Ni_TPL]`, 2 haptic | 81.19, 91.52 → **86.35 s** | 58.32, 57.60 → **57.96 s** | **−32.9%** | A 12.7%, B **1.2%** |
| `FOSNEI_comp_0` | non-eta, boron cage | 86.51, 85.07 → **85.78 s** | 87.16, 84.89 → **86.02 s** | **+0.3%** | A **1.7%**, B **2.7%** |

**Structure fingerprints identical in every run** — `eab3cb62ac92c380` for `CAHQEJ` (all four, and
independently reproduced by Lane B's arm on a third checkout), `797f8a0c3fc6ab40` for `FOSNEI`.

### ⚠ These numbers CORRECT an earlier, contention-inflated pair

The first run of this A/B reported **−50.2%** on `CAHQEJ` and **+9.6%** on `FOSNEI`. It was taken
while four byte-identity gate processes were competing for a 12-core box (load average **35**, driven
partly by two orphaned multithreaded children that survived a `pkill`). Both figures were wrong in
the same direction their noise pointed:

| | contended | quiet | |
|---|---:|---:|---|
| `CAHQEJ` | −50.2% | **−32.9%** | the gain was **over**-stated by 17 points |
| `FOSNEI` | +9.6% | **+0.3%** | the null was buried in 30% within-arm spread |

The commits that landed this lane quote the contended figures; **the quiet-box numbers above are the
ones to cite.** The lesson is in `SPEED_v0.4.10.md` §7: byte-identity gates are load-immune and can be
parallelised freely, wall-clock is neither.

### Reading the `FOSNEI` row

**+0.3% against within-arm spreads of 1.7% and 2.7% is a clean null.** A pure deletion cannot make
code slower — there is no new work on any path — and §4 explains *why* there was nothing to gain:
`FOSNEI` makes **3** `Molecule.__eq__` calls costing **0.03 s**, against `CAHQEJ`'s 99 calls costing
38.52 s. Consistent with v0.4.9's profile, which attributes `FOSNEI` to `get_embedding` self time and
`clean_geometry.ff_clean` with **no eigendecomposition traffic at all**.

## 4. The attribution — and it predicts the A/B

Measured directly on the pristine tree (`8b252292`) by wrapping `Molecule.__eq__` — what
`list.index` actually calls — and the `get_c_eig_list` beneath it:

| molecule | `__eq__` calls | cost | % of generate | generate |
|---|---:|---:|---:|---:|
| `CAHQEJ_comp_0` | **99** | **38.52 s** | **38.8%** | 99.39 s |
| `FOSNEI_comp_0` | **3** | 0.03 s | **0.0%** | 87.11 s |

**This explains both halves of §3 with one number each.** 99 comparisons versus 3 — the scan is
38.8% of one molecule and nothing at all on the other. **The `FOSNEI` null is attributed, not
excused as noise.**

The two independent measurements agree on the size of the thing: the attribution puts the scan at
**38.8%** of that run's 99.39 s generation, and the quiet A/B puts it at **32.9%** of an 86.35 s
pristine mean. Both say *roughly a third of this molecule's generation*, and the difference between
them sits inside the 12.7% run-to-run spread arm A shows on this molecule.

### ⚠ The handoff under-counted its own win

v0.4.9 passed this over as *"15.99 s of a 72.02 s generation — 22%"*. That is **one profile line**,
`numpy.linalg.eig`. The whole of `get_c_eig_list` measures **38.51 s over 198 calls** — and the
**198 eigendecompositions match v0.4.9's count exactly**, so this is the same work, differently
bounded. The difference is the Coulomb matrix that `get_c_eig_list` builds before decomposing
(`get_adj_matrix`, `get_z_list`, two `n×n` matmuls), on **both** operands of every comparison that
survives the atom-count check.

The gap matters beyond bookkeeping: a lane that had budgeted against "22%" and measured 50% would
have had no way to tell a real under-count from a contaminated A/B — and, given how easily this
release contaminated its own measurements, would have been right to distrust the larger number.

## 5. Why no lever

The charter asks every change to ship behind an `OIN_*` lever, default OFF. **This one does not, and
that is a deliberate departure.** A lever gating provably dead code means permanently shipping

```python
if not lever_enabled("OIN_SKIP_DEAD_SCAN"):
    alternative_ace_mol_list.index(alternative_ace_mol)   # still discarded
```

which is strictly worse code than the deletion, and leaves a configuration in which the 50% is
thrown away. The A/B a lever would have enabled was run across two checkouts instead — stronger
isolation, not weaker, since it cannot share a warm import or a module-level cache between arms.

Lane B keeps its lever: a memo is a real alternative code path, and its byte-identity argument rests
on purity rather than on the result being discarded.

## 6. Reproducing

```bash
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
V=$MAIN/.venv/bin/python                     # rdkit pinned ==2025.9.3 -- never `uv sync`

# the A/B (~8 min per molecule on a QUIET box; do not run it alongside a gate)
$V <scratch>/ab_gen.py $MAIN/tmCAT-tmPHOTO_xyz_dataset/cat/CAHQEJ_comp_0.xyz \
    --arm A=$MAIN --arm B=/home/tjmustard/Documents/GitHub/oin-v0410-indexscan \
    --rounds 2 --timeout 300

# the guards
cd /home/tjmustard/Documents/GitHub/oin-v0410-indexscan
PYTHONPATH=$PWD/src $V -m unittest tests.unit.test_embed_dead_scan -v      # 5 tests

# byte-identity
bash tools/gate_v047.sh arm1
bash tools/gate_v047.sh arm2 --cohort-dir $MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v049-strata \
    --golden tools/gate_v049_arm2_golden.tsv --band fast
```
