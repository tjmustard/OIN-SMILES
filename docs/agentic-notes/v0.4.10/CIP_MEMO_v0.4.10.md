# v0.4.10 · Lane B — memoising the CIP re-parse, and closing v0.4.9's ε

> **`VAFMIA_comp_0`: 81.89 s → 10.87 s. −86.7%. Byte-identical.**
>
> That molecule is not an arbitrary example. It is the one that set v0.4.9's budget
> ε — the reason the enforced bound holds to **2.09×, not 1.0×** — because the deadline
> can decline to *start* the next `accept_fn` call but cannot interrupt the one running,
> and the one running was ~24 s of `_reparse_cip_label_once`.

## 1. What was slow

`chirality._reparse_cip_label_once(smiles, probe, fill_deficit)` — one `rdCIPLabeler` pass over a
SMILES, optionally H-filling open valences first. v0.4.9 profiled it at **2.43 s per call, 32 calls,
77.78 s of self time = 99% of `VAFMIA_comp_0`'s generation**
(`docs/agentic-notes/v0.4.9/BUDGET_BOUND_v0.4.9.md` §2).

It is reached two ways, and **both matter**:

| path | side | who calls it per conformer |
|---|---|---|
| `accept_fn` → `_reencode_key_matches` → `build_contract_mol` → `_template_sp3_label` | generator | ARM 2 |
| `chirality.recover()` → `_reparse_aromatic_cip_label` | **encoder** | ARM 1 |

The encoder path is easy to miss — the charter frames this release as generator-side work — and a
memo that was transparent on one side and not the other would sail through a single-arm gate. Both
arms are run below.

## 2. Why a memo is byte-identical, in the form this release requires

> *The removed computation was a pure function of inputs already available, so the output is
> identical and only the wall-clock changes.*

All three arguments are immutable scalars (`str`, `int`, `bool`). The body builds a **fresh**
`Chem.Mol` from the SMILES on every call, touches no module or global state, and returns an
interned `str` or `None`. A hit therefore returns precisely what a miss would compute. Same
sentence as `compute_chg_and_bo_pulp`'s topology memo, which is this project's precedent for the
argument.

**Checked, not assumed** — `tests/unit/test_cip_reparse_memo.py`, 11 tests, pins:

- warm equals cold on every branch of the worker (plain sp3 both chiralities, aromatic-adjacent,
  fused ring, metal-stripped open-valence donor, and two inputs that legitimately return `None`);
- **the key is complete** — `O[C@:99]([O])(C)CC` returns `None` with `fill_deficit=True` and `"R"`
  with `False`, run in *both* orders on a cleared cache, so dropping that key component fails here
  rather than in a corpus sweep. Same for `probe` and for enantiomer collision on `smiles`;
- lever default OFF produces **zero cache traffic**, and `OIN_MEMO_CIP_REPARSE=0` *disables* — the
  `"0"`-is-truthy trap that cost this project 23 test failures across two promotions;
- the real caller `_reparse_aromatic_cip_label` agrees across the lever.

## 3. Why it hits rather than merely being safe

`accept_fn` runs **per conformer**. The key is a SMILES derived from the OIN *template* or from a
metal-free fragment — **neither carries coordinates**. So every conformer of a molecule asks the
same question, and the repeat traffic is structural rather than incidental. That is the whole
mechanism, and it is why the saving is large on a molecule with many conformers and an expensive
CIP centre.

## 4. Measured

`ab_gen.py`, two rounds, arms **alternated** within each round so a drifting box shows up as
within-arm spread rather than as a fake between-arm effect. Same checkout for both arms; only the
lever differs. One fresh interpreter per (molecule, arm).

### `VAFMIA_comp_0` — `[Cu_LIN]`, bis-adamantyl NHC

| arm | runs | mean | spread |
|---|---|---:|---:|
| OFF | 83.55, 80.24 | **81.89 s** | 4.0% |
| ON | 11.03, 10.71 | **10.87 s** | 2.9% |

**−71.02 s = −86.7%**, against a **0.28%** noise floor — ~310× the floor.
**Structure fingerprint `99f1650d65f245f1` in all four runs.**

> **It crosses the 30 s line.** VAFMIA spent 62.8 s against a 30 s budget in v0.4.9's A/B — the
> single worst overrun and the definition of that release's ε. At 10.87 s it is now comfortably
> inside. That is one molecule, not the goal: v0.4.9's tail had **11** molecules over budget in
> both arms, and this closes the one that was worst, not all eleven.

## 5. What this does NOT claim

- **Not a corpus number.** The cost is **bimodal by molecule** — v0.4.9 measured the same
  `accept_fn` re-encode at **62% of VAFMIA and 0.8% of `FOSNEI_comp_0`**. A lever aimed at whichever
  function profiled expensive last optimises one molecule class and measures nothing on the corpus;
  v0.4.9 shipped that mistake once and recorded it. The stratified tail sample in
  `SPEED_v0.4.10.md` is the honest multi-molecule figure.
- **Not `max(elapsed_s) < 30 s`.** See above.
- **Not an encode-side optimisation lane.** The encoder benefit is a side effect of the same memo,
  not a separate change; the forked-resonance R3 regime that dominates encode cost is v0.4.11.

## 6. What shipped

Behind **`OIN_MEMO_CIP_REPARSE`, default OFF**. Unset, `_reparse_cip_label_once` dispatches
straight to the un-memoised worker and the cache is never touched.

| change | file |
|---|---|
| worker renamed `_reparse_cip_label_once` → `_reparse_cip_label_once_uncached`; new `lru_cache` slot; dispatcher keeps the old public name | `src/oinsmiles/core/chirality.py` |
| `_reparse_cip_memo_clear()` / `_reparse_cip_memo_info()` — isolation + telemetry, mirroring `_ac2bo_memo_clear()` | `src/oinsmiles/core/chirality.py` |
| lever registered with its rationale and the bimodality warning | `src/oinsmiles/oin/levers.py` |
| 11 guards | `tests/unit/test_cip_reparse_memo.py` |

Bounded at `_CIP_REPARSE_MEMO_MAX = 2048` so a long single-interpreter sweep cannot grow it without
limit. Cross-molecule retention is safe by §2's purity argument; the clear hook exists so a
per-molecule gate can guarantee isolation regardless.

## 7. Reproducing

```bash
MAIN=/home/tjmustard/Documents/GitHub/OIN-SMILES
V=$MAIN/.venv/bin/python                     # rdkit pinned ==2025.9.3 -- never `uv sync`
W=/home/tjmustard/Documents/GitHub/oin-v0410-cipmemo

# the A/B  (~6 min, quiet box)
$V <scratch>/ab_gen.py $MAIN/tmCAT-tmPHOTO_xyz_dataset/cat/VAFMIA_comp_0.xyz \
    --arm OFF=$W --arm ON=$W --env ON:OIN_MEMO_CIP_REPARSE=1 --rounds 2 --timeout 300

# the guards
cd $W && PYTHONPATH=$PWD/src $V -m unittest tests.unit.test_cip_reparse_memo -v   # 11 tests

# byte-identity, BOTH arms -- the encoder path is not optional here (§1)
cd $W && OIN_MEMO_CIP_REPARSE=1 bash tools/gate_v047.sh arm1
cd $W && OIN_MEMO_CIP_REPARSE=1 bash tools/gate_v047.sh arm2 \
    --cohort-dir $MAIN/tmCAT-tmPHOTO_xyz_dataset/cohort-v049-strata \
    --golden tools/gate_v049_arm2_golden.tsv --band fast
```
