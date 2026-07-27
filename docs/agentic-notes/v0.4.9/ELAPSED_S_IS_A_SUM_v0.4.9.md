# v0.4.9 · `elapsed_s` is a SUM — this release's own premise, refuted

> **v0.4.9 was chartered on one number: "759.9 s against a 300 s budget. That single number is
> this release's justification: the budget is not a budget."**
>
> That number is arithmetic on a sum. `metrics.elapsed_s` accumulates **up to three separately
> supervised attempts**, each under its own SIGKILL. Measured over the 5000-molecule sweep:
> **every single-attempt run finished within 0.2 s of the 300 s budget.** The harness's kill
> enforces to ε ≈ 0.2 s. 759.9 s is 300 + 300 + 160.
>
> The release proceeds. Its justification does not.

## 1. Where the sum comes from

`tools/test_dataset_roundtrip.py` runs each molecule through up to three supervised attempts:

| stage | line | budget |
|---|---|---|
| PASS 1 `UFF_1` | `:843` | `args.mol_timeout`, SIGKILL via `_supervise()` |
| PASS 2 tier 1 (`ensemble_size=1`) | `:914` | `args.mol_timeout`, again |
| PASS 2 tier 5 (`ensemble_size=5`) | `:948` | `args.mol_timeout`, again |

and **adds their wall-clock into one field** (`:850`, then `elapsed_s +=` at `:916` and `:951`).
The ceiling is therefore `3 × mol_timeout = 900 s`, and 759.9 s sits comfortably inside it.

## 2. The measurement that settles it

Split the 5k sweep by `tier_passed` — which records *how many* attempts a molecule consumed:

| `tier_passed` | n | median | max | `> 310 s` |
|---|---:|---:|---:|---:|
| `UFF_1` — **one** attempt | 4658 | 6.72 s | **300.2 s** | **0** |
| `g-xTB_1` — two attempts | 2 | 330.15 s | 332.8 s | 2 |
| `None` — failed, ran all three | 340 | 300.33 s | **759.9 s** | 83 |

**4658 single-attempt runs, worst case 300.2 s against a 300 s budget.** Not one exceeds it by
more than 0.2 s. Every row above 310 s consumed more than one attempt.

```bash
V=$PWD/.venv/bin/python; export PYTHONPATH=$PWD/src
$V - <<'PY'
import json, glob, os, statistics as st
D = "tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/individual_reports"
rows = []
for f in glob.glob(os.path.join(D, "*.json")):
    r = json.load(open(f))
    e = r.get("metrics", {}).get("elapsed_s")          # ⚠ NESTED
    if e is not None:
        rows.append((e, r.get("tier_passed")))
for tier in {t for _, t in rows}:
    sub = sorted(e for e, t in rows if t == tier)
    print(f"{str(tier):10s} n={len(sub):5d} median={st.median(sub):7.2f} max={sub[-1]:7.1f} "
          f">310s={sum(1 for x in sub if x > 310)}")
PY
```

## 3. What survives

**The in-process advisory timeout is real** — but the evidence is the code and two direct-call
probes, never the 759.9 s figure.

`generator3d/__init__.py:350` computes the deadline; `:457` and `:529` check it **at the top of
the attempt loop only**. Everything inside one iteration runs unbounded:

- the `accept_fn` re-encode — a full `XYZToSMILES().convert()`, **measured 48–57 s per call** on
  an eta/haptic conformer (`generation/metallogen_adapter.py:2118`), invoked *after* the deadline
  check;
- the PuLP/CBC bond-order solve in `alt_cache` priming — `prob.solver.actualSolve(prob)` with
  **no `timeLimit`** (`generator3d/utils/compute_chg_and_bo_pulp.py:254-260`), one subprocess per
  solve, counter-attributed as the boron cost sink;
- `AllChem.EmbedMolecule`;
- `_select_by_geometry` (`metallogen_adapter.py:2204`), which runs *after* the loop returns.

Two independent probes measured the consequence directly, with no SIGKILL in the way: 60 s asked
→ 60.7–137.9 s (eta sample, `GOHWOQ` **2.3×**) and 60.0–172.8 s (boron set, **2.9×**).

**So: the library's contract is advisory. The harness's kill is not.** A consumer that is not the
harness has no bound at all — which is the honest statement of the problem v0.4.9 exists to fix.

## 4. The prize, correctly sized — 13× the charter's estimate

The charter justified the bound on the boron cap-burners' **~2.1 CPU-h per sweep**. That is 8% of
the real waste.

| band | n | CPU-h | % of sweep | honest passes | CPU per pass |
|---|---:|---:|---:|---:|---:|
| `< 30 s` | 4006 | 8.67 | 15.8% | 3380 (84.4%) | 0.2 min |
| `30–60 s` | 292 | 3.45 | 6.3% | 136 (46.6%) | 1.5 min |
| `60–300 s` | 411 | 15.36 | 28.0% | 112 (27.3%) | 8.2 min |
| **`≥ 300 s`** | **291** | **27.31** | **49.8%** | **3 (1.0%)** | **546 min** |
| total | 5000 | **54.8** | | 3631 | |

**Half the sweep's compute buys three molecules.**

## 5. And the bound is cheaper than anyone assumed

**93.1% of all honest passes already finish in under 30 s.** So capping per-molecule wall-clock
costs a known, small, one-time number of passes:

| per-molecule cap | CPU-h recovered / sweep | passes lost | `byte_exact` cost |
|---:|---:|---:|---:|
| 300 s | 3.06 | 3 | 0.06 pts |
| 120 s | 21.33 | 52 | 1.04 pts |
| 60 s | 30.97 | 115 | 2.30 pts |
| **30 s** | **37.84** | **251** | **5.02 pts** |

The roadmap treated Goal A (100% `byte_exact`) and Goal B (`max(elapsed_s) < 30 s`) as one
problem because "78.8% of failures never test the notation". That is true of *failures*; it is
not true of *passes*. **The goals are separable, at a price of about five points**, and a 5k
sweep under a 30 s bound costs **~17 CPU-h instead of 55** — which makes live sweeping an
affordable instrument again, rather than something v0.4.8 had to work around with an offline
re-score.

## 6. Two things the charter did not know

**`--mol-timeout` never reaches the generator.** `tools/test_dataset_roundtrip.py:808` hardcodes
`timeout_val = 30 if args.quick else 300` and passes *that* to `OIN3DGenerator(timeout=)`;
`--mol-timeout` feeds only `_supervise()`'s SIGKILL. **There is no way to ask the harness for a
different generator budget**, so the A/B this release needs is not currently expressible. Fixing
the plumbing is a prerequisite, not scope creep.

**15 molecules burn the full budget inside the ENCODER.** Every one reads
`TimeoutException at UFF_1: exceeded 300s while encoding (hard kill)` — 1.25 CPU-h, and exactly
the `encode_fail` bucket (15, 0.30%). **No generator-side bound can reach them**: they never get
as far as the generator. Only a budget enforced at or above `OIN3DGenerator.generate()` — i.e.
the per-molecule form — covers them. That is an argument for the per-molecule design, and a limit
that must be stated wherever a generator-only bound is claimed to deliver Goal B.

## 7. Why this is recorded rather than quietly fixed

Four of this project's releases have ended by refuting their own plan. The failure mode this one
illustrates is specific and worth naming:

> **A derived aggregate was read as a direct observation.** `elapsed_s` looks like a duration. It
> is a sum over retries, and nothing in its name says so. The same field had already produced one
> wrong claim by being read from the top level instead of from `metrics` — that trap is documented
> everywhere in this repo. This is the second trap in the same field, and it survived the first
> audit because everyone was busy remembering the first one.

**Both traps, stated together:** `metrics.elapsed_s` is **nested**, and it is a **sum across up
to three supervised attempts**. A runtime claim must say which it means.
