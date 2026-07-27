# v0.4.10 · Lane A — deleting one line that eigendecomposed a matrix and threw it away

> **`CAHQEJ_comp_0`: 118.44 s → 59.00 s. −50.2%. Byte-identical.**
> **`FOSNEI_comp_0`: no measurable effect.**
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
| `CAHQEJ_comp_0` | eta, `[Ni_TPL]`, 2 haptic | 117.94, 118.95 → **118.44 s** | 61.47, 56.53 → **59.00 s** | **−50.2%** | A 0.9%, B 8.4% |
| `FOSNEI_comp_0` | non-eta, boron cage | 83.86, 94.36 → **89.11 s** | 84.77, 110.63 → **97.70 s** | +9.6% | A 12%, **B 30%** |

**Structure fingerprints identical in every run** — `eab3cb62ac92c380` for `CAHQEJ` (all four, and
independently reproduced by Lane B's arm on a third checkout), `797f8a0c3fc6ab40` for `FOSNEI`.

### Reading the `FOSNEI` row honestly

**+9.6% is noise, not a slowdown.** A pure deletion cannot make code slower; there is no new work on
any path. The within-arm spreads (12% and 30%) are larger than the between-arm difference, which is
the definition of an unresolved measurement — these two rounds ran while four gate processes were
competing for the box. The defensible statement is **"no measurable effect on `FOSNEI`"**, and it is
consistent with v0.4.9's profile of that molecule, which attributes its cost to `get_embedding` self
time and `clean_geometry.ff_clean` and shows **no eigendecomposition traffic at all**.

The `CAHQEJ` row is the opposite: within-arm spread of 0.9% on the arm that matters, against a
between-arm difference of 50%.

## 4. ⚠ The handoff under-counted its own win — 22% predicted, 50% measured

v0.4.9 handed this over as *"15.99 s of a 72.02 s generation — 22%"*. That figure is **one profile
line**: `numpy.linalg.eig` via `chem.get_c_eig_list`. The scan's true cost is everything inside
`Molecule.__eq__`, which also includes building the Coulomb matrix (`get_adj_matrix`, `get_z_list`,
two `n×n` matmuls) on **both** operands of every comparison that survives the atom-count check.

The gap matters beyond bookkeeping: a lane that had budgeted its effort against "22%" and measured
50% would have had no way to tell a real under-count from a contaminated A/B, and would have been
right to distrust the larger number. Attribution is in `SPEED_v0.4.10.md`.

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
