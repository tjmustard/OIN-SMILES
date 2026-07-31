# The remaining 22.84 points, partitioned by **what has actually been measured to reach them**

**Derived 2026-07-30 from the v0.4.14 baseline sweep (n=5000, `byte_exact` 77.16%) plus v0.4.15's
six A/B arms.** This is the input to the v0.4.16–v0.4.18 charters. Every row is a measurement or
is explicitly marked UNMEASURED — no row is a hypothesis wearing a bucket name.

## The map

| block | n | pts | selection-reachable? | evidence |
|---|---:|---:|---|---|
| `rdkit_canonical` | 113 | 2.26 | 🟢 **38.1% (43 mol, 0.86 pts)** | v0.4.15 L2 arm |
| `slot_renumber` (non-enantiomer part) | 51 | 1.02 | 🟢 partial (5 mol) | v0.4.15 L2 arm |
| `slot_renumber` → **enantiomers** | 201 | 4.02 | 🔴 **NO — 1 of 201 (0.5%)** | v0.4.15 L2 arm + telemetry |
| `structural` → `DETACHED` | 301 | 6.02 | 🔴 **NO — 0 of 289** | v0.4.15 L1 arm |
| `structural` → `INTACT` (MEDZUR) | 110 | 2.20 | ⚪ **UNMEASURED** | — |
| `structural` → `BOUNDARY` | 62 | 1.24 | ⚪ **UNMEASURED** | — |
| `structural` → `NO_STRUCTURE` | 11 | 0.22 | 🔴 nothing generated | audit |
| `hard_fail` | 266 | 5.32 | 🔴 **262 produce NOTHING** | audit |
| `facmer_divergent` | 15 | 0.30 | ⚪ unmeasured (11 also `DETACHED`) | audit |
| `encode_fail` | 12 | 0.24 | 🔴 encoder floor | — |

**Reachable today: 0.96 pts** (built, measured, sitting behind `OIN_ACCEPT_STRING_EXACT` at 4.00×
runtime).
**Proven NOT reachable by selection: 10.04 pts** (the 201 enantiomers + the 301 `DETACHED`).
**UNMEASURED: 3.74 pts** (MEDZUR 110 + BOUNDARY 62 + `facmer` 15).

## 🔴 Correction to v0.4.14: `rdkit_canonical` is the most reachable block on the board

v0.4.14 re-filed `rdkit_canonical` as **80.7% η-set denticity drift** and wrote it out of the
ladder:

> **`rdkit_canonical` (114 / 2.28 pts)** — re-filed as η-set denticity drift. Not reachable by
> canonicalization; **do not re-propose a string fix.**

**v0.4.15's Lane 2 recovered 43 of its 113 — 38.1%, the highest rate of any block measured.**

The two statements are compatible and the distinction is the point:

* v0.4.14's claim was about the **encoder**. It stands: no canonicalization of the emitted string
  fixes an η-set that genuinely differs.
* What is new is that it is reachable from the **generator**. The molecules were being handed a
  conformer whose η-set did not match; filling the pool further finds one that does. **A different
  conformer, not a different string.**

⚠ **The generalisable error is the blanket word "reachable".** "Not reachable by canonicalization"
was recorded, and then read for a release as "not reachable". A block's reachability is a property
of *a mechanism*, not of the block — and this project has now made that mistake in both directions
in consecutive releases (v0.4.14 wrote off a reachable block; v0.4.15's charter aimed two lanes at
unreachable ones).

## The other correction: `slot_renumber` is mostly not reachable either

| `key_equal` subclass | n | reachable now | rate |
|---|---:|---:|---:|
| `rdkit_canonical` | 113 | **43** | **38.1%** |
| `slot_renumber` | 252 | 5 | **2.0%** |

`slot_renumber` is 5.04 points and **201 of its 252 are the enantiomer class**, which is
construction-blocked. So the headline "`slot_renumber` is the largest single encoder block" —
true since v0.4.11 — describes a block that is 80% unreachable by anything the encoder or the
selector can do.

## What this means for the ladder

1. **The gap is now majority-CONSTRUCTION.** 10.04 of 22.84 points are proven unreachable by
   selection, and another 5.32 produce no structure at all. **15.36 points (67%) require the
   generator to build something different**, not to choose differently or to emit differently.
2. **0.96 points are sitting behind a lever that is built and measured**, blocked only on the
   runtime trade.
3. **3.74 points have never been characterised** — MEDZUR (110), BOUNDARY (62), `facmer` (15).
   That is now the largest *unknown*, and the project's own rule applies: *characterize a block's
   members before pointing a lane at it.*

## Reproduce

```bash
V=.venv/bin/python; export PYTHONPATH=$PWD/src
D=tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep
$V tools/roundtrip_bucket_report.py --results-dir "$D" --score honest   # 3858 = 77.16%
$V tools/attach_class_audit.py --results-dir "$D"                       # structural split
# the reachability columns come from measurements/v0.4.15/lane{1,2}_pop_*.json
```
