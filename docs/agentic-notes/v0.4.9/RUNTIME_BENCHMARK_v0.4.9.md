# v0.4.9 · the frozen runtime benchmark — provenance, strata, and blind spots

> **What Goal B (`max(elapsed_s) < 30 s`) is measured on, from here on.**
>
> 328 molecules, stratified across four runtime bands × eta/non-eta, drawn from **one** source
> by a deterministic rule, with every row labelled by its honest round-trip class.
> `tools/gate_v049_arm2_golden.tsv` freezes 325 of them.

## 1. Why the inherited cohort was replaced

`tools/gate_v047_arm2_golden.tsv` is the top 100 molecules by `elapsed_s` that were byte-exact
in `results-v0.4.5-sweep-partial-2697mols`. Its author flagged it provisional. Three
independent problems, only the first of which was known:

**(a) Wrong shape for a p100 goal.** Cut at 93.06 s, it cannot see the 30–93 s band — where
most of the 994 over-30 s molecules live — and it has **no fast control at all**, so a change
that speeds the tail while slowing everything else measures as a pure win.

**(b) 🔴 59% of it is a false positive.** Cross-tabulated against v0.4.8's honest baseline:

| honest class | n |
|---|---:|
| `byte->byte` | 28 |
| **`byte->FAIL`** | **59** |
| `byte->key` | 4 |
| `no_structure` | 8 |
| `key->FAIL` | 1 |

It selected on `smiles_1 == smiles_2` — the **scored** verdict, which reads the generator's own
bond graph. v0.4.8 showed that over-states `byte_exact` by 10.34 points corpus-wide and by far
more among slow haptic molecules, which is exactly what a slow cohort is made of. Eight members
additionally produced no structure at all in the 5k sweep.

**(c) It was drawn from a different, smaller sweep** than the one every current number comes
from.

None of this makes the v0.4.7 gate *wrong* — a byte-identity gate detects change, and a frozen
false-positive string still detects change. It makes it **unrepresentative** of the corpus whose
p100 is the goal, and dangerous to read as a correctness statement.

## 2. The draw

| | |
|---|---|
| **source** (ONE dir, never intersected) | `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest` |
| **why that one** | it is `results-v0.4.6-sweep` plus the honest fields — the only dir carrying both `metrics.elapsed_s` and `honest_class` |
| **source sweep commit** | v0.4.6-era generation, honest-rescored at v0.4.8 |
| **selection predicate** | **none on pass/fail.** Only a usable `metrics.elapsed_s`. |
| **sampling** | systematic on the elapsed-sorted rank within each cell, `round(i·(n−1)/(quota−1))` — deterministic, no seed |
| **tool** | `tools/select_runtime_strata.py` |
| **date** | 2026-07-27 |
| **manifest** | `spec/handoffs/v0.4.9/cohort_v049_strata.json` (gitignored; the golden below is the durable artifact) |
| **cohort dir** | `tmCAT-tmPHOTO_xyz_dataset/cohort-v049-strata` (328 symlinks) |

**No pass/fail predicate is applied, deliberately.** A runtime benchmark must contain the
molecules that fail — they are the slow ones. Requiring `status=success` is how the 30–300 s
band went missing the first time. Every row instead carries its `honest_class`.

### Strata

| band | corpus n | selected eta | selected non-eta | total |
|---|---:|---:|---:|---:|
| `fast` (< 30 s) — **the control** | 4006 | 30 | 60 | 90 |
| `b30_60` (30–60 s) — the band v0.4.7 omitted | 292 | 40 | 40 | 80 |
| `b60_300` (60–300 s) | 411 | 45 | 40 | 85 |
| `b300_plus` (≥ 300 s) — budget already violated | 291 | 35 | 35 | 70 |
| `unbaselined` pins | — | — | — | 3 |
| | 5000 | 150 | 175 | **328** |

Quotas are **not** proportional to corpus share. The `fast` cell is a control — it exists to
catch "the tail got faster and everything else got slower", which nothing in the v0.4.7 cohort
could detect — while the three slow cells are what Goal B is about and need enough members to
resolve a shift.

### Composition, by honest class

| honest class | n | | | n |
|---|---:|---|---|---:|
| `byte->byte` | 128 | | `key->FAIL` | 13 |
| `no_structure` | 71 | | `fail->fail` | 2 |
| `byte->FAIL` | 71 | | `FAIL->byte` | 1 |
| `key->key` | 21 | | `FAIL->key` | 1 |
| `byte->key` | 17 | | *(unbaselined)* | 3 |

249 `success`, 76 `failed`, 3 unbaselined. Compare the v0.4.7 cohort's 28/59 — this one is not
dominated by false positives, and it deliberately contains the failures.

### Pins

Five molecules are included regardless of what the stride picks, because a benchmark that loses
its own worst case to an arithmetic accident is not a benchmark.

| pin | why | status |
|---|---|---|
| `FOSNEI_comp_0` | worst observed `elapsed_s` (759.9 s) | in source, frozen |
| `RAWJEG_comp_0` | boron cage that **does** assemble (~2.5 s) — refutes the blanket fast-fail | in source, frozen |
| `GOHWOQ_comp_0` | the 2.3× advisory overrun (60 s asked, 137.9 s spent) | **not in the source sweep** |
| `ULODUU_comp_0` | second boron success (TET); budget-dependent at ~61.8 s | **not in the source sweep** |
| `XIQKOY_comp_0` | boron two-point proof (OFF fails 0.87 s, ON runs past 340 s) | **not in the source sweep** |

Three of the release's own named fixtures **are not in the frozen seed-42 5k draw at all.**
They are in the cohort dir and are measured; they have no frozen golden row, and
`--allow-unbaselined` reports that rather than inventing one.

## 3. What each golden row actually gates

325 rows. **A hash is frozen only where byte-identity is a stable property of the code.**

| column-3 value | n | gated on | why |
|---|---:|---|---|
| a real `sha256` | 271 | `sha1` **and** `sha2` | ordinary byte-identity |
| `NO_STRUCTURE@300s` | 45 | `sha1` only | the source run **timed out**. Whether it generates now is a fact about the box — `ULODUU` assembles at a 60 s cap and not at 30 s, the very reason the boron fast-fail was refuted. Freezing an empty hash would make every future speedup a MISMATCH. |
| `NO_STRUCTURE_DET` | 5 | `sha1` **and** "still produces nothing" | the source run **raised** (e.g. `GURKUA` in 1.2 s on `UncoordinatedFragmentError`). That IS code-determined, so it is gated strictly: producing a structure now is a MISMATCH. |
| `NO_ENCODE@300s` | 4 | nothing — observation only | the **encoder** was hard-killed at the budget. No string exists in either column. |

The verdict line reports the four counts separately, and the gate **refuses to pass** if zero
molecules were actually gated — a green verdict over an all-sentinel selection certifies nothing.

Columns 7–8 carry `band` and `honest_class` and are never compared. They exist so nobody reads a
green ARM 2 as "the chemistry is right": **71 of these rows are `byte->FAIL`**, and the gate is
deliberately blind to that.

## 4. Reproducibility, and the noise floor

Two consecutive `--band fast` runs (90 molecules) on an unchanged tree:

| run | verdict | wall-clock |
|---|---|---:|
| 1 | PASS — 87 compared + 3 deterministic-no-structure | **277.01 s** |
| 2 | PASS — 87 compared + 3 deterministic-no-structure | **277.79 s** |
| | | **Δ 0.78 s = 0.28%** |

Byte-level: the per-molecule TSV rows are **identical** between runs (`diff` on columns 1–6 is
empty), so the verdict is not merely stable, it is the same bytes.

> ### The noise floor is **0.28%**
>
> Every later runtime claim must clear it. Without a measured floor, a 10% "improvement" is
> unfalsifiable — and this project has already published a runtime figure that turned out to be
> an artifact.

Both runs were back-to-back on the same box under comparable load. **A figure taken on a
contended box is a ratio within a run, never an absolute** — 0.28% is the floor for
*back-to-back on a quiet box*, and nothing else. Do not compare an idle-box arm against one
measured during a sweep; two of this project's three published overrun ranges were taken
alongside other work, which is exactly why neither is a clean distribution.

## 5. 🔴 Blind spots — what this benchmark cannot see

Naming these is part of the deliverable.

1. **A narrow regression in fast molecules.** The `fast` cell samples 90 of 4006 (**2.2%**). A
   change that slows only one metal, geometry, or ligand class inside the fast band will very
   likely be missed. The control detects *systematic* slowdown, not *targeted* slowdown.
2. **Anything outside the frozen 5k draw.** The corpus is 25,197 unique basenames; the source
   sweep is a 5000-molecule seed-42 sample. Three of this release's own pins are already proof
   that named molecules fall outside it.
3. **Correctness.** ARM 2 gates strings, never chemistry. 71 rows are known false positives and
   pass the gate by construction. Use `tools/roundtrip_bucket_report.py --score honest` for
   correctness; never this.
4. **Budget-dependent generation, by design.** The 45 `NO_STRUCTURE@300s` rows will not fail the
   gate if they start or stop generating. That is deliberate (§3) and it means **a change that
   flips molecules across the budget boundary is invisible here** — it must be measured as a
   count, not gated.
5. **`elapsed_s` bands are per-molecule TOTALS**, summed over up to three tiers. A molecule can
   move between bands because its retry count changed rather than because anything got faster.

## 6. Two tooling defects found while building this

Both were latent in the shipped v0.4.7 gate.

**The gate picked an arbitrary interpreter.** Its fallback globbed
`$(dirname $REPO)/*/.venv/bin/python` because worktrees have no `.venv`. Run from a v0.4.9
worktree under `~/Documents/GitHub` that matched **`EtaCatalysis/.venv`** — an unrelated project
with no rdkit — and then **`EtaTMCSMILES/.venv`**, which has rdkit **2025.09.2** against this
project's pinned **2025.9.3**. The first fails loudly but is misdiagnosed (a short `#DONE` reads
as "the run was truncated"); the second is worse, because **a byte-identity gate on a different
rdkit reports MISMATCHes that look like code regressions.** Now resolved deterministically via
`git rev-parse --git-common-dir`, with the rdkit version checked against the `pyproject.toml`
pin and a hard refusal on drift.

**The gate itself was unbounded.** `--timeout` is only `OIN3DGenerator`'s *advisory* budget, so
ARM 2 had no enforced cap — the same defect Lane 1 exists to fix, inside the tool measuring it.
Added `--hard-timeout` (real `SIGKILL`, default 1.5× the generator budget) with a synthesized
`HARD_TIMEOUT@Ns` row, so a killed molecule is a *result* rather than a missing line that makes
`#DONE` short.

Also added: `--shard i:n` (**1-based**, matching `test_dataset_roundtrip.py`, where `--shard 0:6`
exits 2 and launching 0..5 silently drops a sixth of the cohort) and `--band NAME`, because the
full cohort is ~10 CPU-h serially and nobody runs that twice for a reproducibility check.

## 7. Running it

```bash
cd <checkout>
D=$PWD/tmCAT-tmPHOTO_xyz_dataset

# routine — the fast control, ~4.5 min
bash tools/gate_v047.sh arm2 --cohort-dir $D/cohort-v049-strata \
    --golden tools/gate_v049_arm2_golden.tsv --band fast

# release — the whole cohort, sharded (1-BASED)
for i in 1 2 3 4 5 6; do
  bash tools/gate_v047.sh arm2 --cohort-dir $D/cohort-v049-strata \
      --golden tools/gate_v049_arm2_golden.tsv --shard $i:6 --timeout 300 \
      --out /tmp/v049_arm2_$i.tsv &
done; wait

# re-derive the cohort from scratch (deterministic -- same 328 names)
PYTHONPATH=src .venv/bin/python tools/select_runtime_strata.py \
    --results-dir $D/results-v0.4.8-honest \
    --out spec/handoffs/v0.4.9/cohort_v049_strata.json \
    --names-out spec/handoffs/v0.4.9/cohort_v049_names.txt
PYTHONPATH=src .venv/bin/python tools/build_sweep_cohort.py \
    --names-file spec/handoffs/v0.4.9/cohort_v049_names.txt \
    --dataset-dir $D --out $D/cohort-v049-strata

# rebuild the golden
PYTHONPATH=src .venv/bin/python tools/gate_v047_build_arm2_golden.py \
    --names-file spec/handoffs/v0.4.9/cohort_v049_names.txt \
    --source-a $D/results-v0.4.8-honest \
    --predicate any_encoded --source-budget-s 300 \
    --strata-file spec/handoffs/v0.4.9/cohort_v049_strata.json \
    --allow-unbaselined --out tools/gate_v049_arm2_golden.tsv
```

**ARM 1 is untouched** — 61 fixtures, same golden, same verdict.

The v0.4.7 golden also still rebuilds **bit-identically** under the default
`--predicate byte_exact_scored`; that is checked, not assumed.
