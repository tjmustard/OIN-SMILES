# v0.4.16 Lane 1 — `OIN_ACCEPT_STRING_EXACT`, priced · **HOLD: the bound does not create a knee**

**Verdict: HOLD, and the negative is the finding.** Bounding the search does not buy a favourable
trade — **cost and gain rise together**, so there is no bound where most of the +48 is retained
cheaply. The lever stays default-OFF, `OIN_STRING_EXACT_BOUND` ships default-unset, and the
release is FLAT at **77.16%**.

## The question

v0.4.15 landed the lever measured at **+48 molecules / +0.96 pts, zero losses over 565** — and at
**4.00× runtime**, `> 30 s` going 30 → 122 on its 365-molecule population (~678 → ~770
corpus-wide). The roadmap targets `byte_exact` 100% **and** `max(elapsed_s) < 30 s`, so that is
close to a wash between the two halves. The cost has a known mechanism: the lever declines to
*stop* the pool, so the pool fills to budget. Re-read from the frozen arm, **the 317 molecules that
never gain consume 16149 s of the lever's 17191 s — 93.9% of the bill buys nothing.**

So: bound the extra search, and find the knee.

## The curve — derived from ONE run, then confirmed

n = 365, all rows scored. **2 excluded** (no incumbent recorded — for those, truncation is not
answer-neutral, so they are unmeasurable by this method and are *not* counted as zeros).

| bound | recovered | of ceiling | total s | `> 30 s` |
|---:|---:|---:|---:|---:|
| 0 | 0 | 0.0% | 4534 | 31 |
| 1 | 11 | 22.9% | 5938 | 34 |
| 2 | 15 | 31.2% | 7528 | 42 |
| **3** | 19 | 39.6% | 8564 | **52** |
| 5 | 23 | 47.9% | 10964 | 69 |
| 9 | 27 | 56.2% | 14516 | 97 |
| 10 | 31 | 64.6% | 14658 | 101 |
| **12** | **38** | **79.2%** | 14975 | 104 |
| 14 | 40 | 83.3% | 15231 | 104 |
| 19 | 42 | 87.5% | 15789 | 109 |
| 25 | 48 | 100.0% | 16212 | 113 |
| ∞ | 48 | 100.0% | 16418 | 117 |

## 🔴 The verdict: the two conditions are mutually unreachable

The owner-set bar: 0 losses **and** ≥ 36 of 48 recovered (≥ 75%) **and** `> 30 s` on the 365 ≤ 52.

| | needs | gives |
|---|---|---|
| **condition 2** (≥ 36 recovered) | bound **≥ 12** | `> 30 s` = **104** ✗ |
| **condition 3** (`> 30 s` ≤ 52) | bound **≤ 3** | recovered **19** ✗ |

**No bound satisfies both.** They are separated by a factor of four in the bound and a factor of
two in each metric.

### Why there is no knee — the shape, not the threshold

A threshold can be argued with. The shape cannot:

```
bounds  0 → 9    +9982 s   for  27 molecules     (the pool's expensive early conformers)
bounds  9 → 12    +459 s   for  11 molecules     ← the ONE efficient segment
bounds 12 → ∞    +1443 s   for  10 molecules
```

To keep **79%** of the gain you pay **89%** of the runtime penalty (104 of 117 extra `> 30 s`
molecules). **The frontier is close to linear**: bounding moves along the v0.4.15 trade rather than
improving it. The charter's hypothesis — that a wasted tail could be reclaimed cheaply — is
**refuted**, and the reason is that the tail is not where the cost is. Each post-incumbent
conformer is a full embed *plus* a full `XYZToSMILES` re-encode, so the first few extra conformers
are the expensive ones.

⚠ Charter candidate 1 — *"verify the pool is not filling past the first string-exact conformer"* —
is answered and it was already true: `early_hit` breaks both fill loops immediately. There was no
wasted tail to reclaim.

## The instrument was validated at three independent points, not asserted

The curve is **derived** from one instrumented run rather than measured per bound. That is only
legitimate because `incumbent_hit` is returned whatever the pool does afterwards. Three checks,
none of which the derivation could fake:

| check | derived | independently measured (v0.4.15 frozen arm) |
|---|---:|---:|
| ceiling (recovered, unbounded) | **48** | **48** — exact |
| bound-0 total runtime | 4534 s | 4294 s (+5.6%) |
| unbounded total runtime | 16418 s | 17191 s (−4.5%) |
| bound-0 `> 30 s` | 31 | 30 |
| unbounded `> 30 s` | 117 | 122 |

And the **live bound-0 wiring gate**, over 40 molecules:

```
gains 0    losses 0    GENERATED output moved 0
```

That last field is the strong form: bound 0 does not merely reach the same *verdict* as lever-OFF,
it produces **byte-identical generated output**. A broken bound and a working one agree at large
*N* and differ at 0 — this is the check that separates them, and it passes.

### 🔴 The live bound-12 arm: 48 of 48, in BOTH directions

The derivation predicts, per molecule, whether bound 12 still recovers it. Run live over **all 48**
hit molecules:

| | derived | live | |
|---|---:|---:|---|
| KEEP (`min_bound` ≤ 12) — must gain | 38 | **38** | agree |
| DROP (`min_bound` > 12) — must **not** gain | 10 | **10** | agree |
| losses | 0 | **0** | |

**48/48, 100%, both directions.** And the 10 dropped molecules read `out_moved=False` — the bound
stopped the search and returned output byte-identical to the OFF arm, which is what a correct
truncation must do.

⚠ **The two-directional part is the whole point, and the first attempt did not have it.** The
population file was sorted by `min_bound`, so when an arm was killed at 30 of 48 every completed
row happened to be a molecule the bound *keeps*. That arm read 30/30 agreement and confirmed
nothing about the bound — an arm containing only cases the mechanism accepts cannot discriminate.
Same shape as v0.4.13's *"ARM 1 is BLIND — 0 of 62 fixtures are fold-movers, so PASS ≠
validation"*. The 18 remaining molecules, 10 of which are the discriminating ones, were re-run
specifically to close it.

Run conditions were verified comparable to the frozen arm before anything was quoted: **0.99×
total, 0.95× median** over the same molecules. Load biases *accuracy* here, not just timing —
the generator timeout is advisory, so starvation shrinks the pool.

## What the owner may want to override

The bar was set **before this data existed**, so it is worth stating the two points it rejects:

- **bound 3** — +19 molecules (**+0.38 pts**) for +21 over-30 s (corpus ~699, inside the ~700
  budget). The cheapest real gain on the board; fails only the ≥ 75% retention rule.
- **bound 12** — +38 molecules (**+0.76 pts**) for +73 over-30 s (corpus ~751). The best accuracy
  available under any bound; fails only the runtime rule.

Both are defensible if the two halves of the target are weighted differently. **Neither is
recommended here**, because the measured shape says bounding is not the lever that resolves this
trade — and shipping a default that trades ~1 point of accuracy against ~1.5 of runtime should be
a deliberate act, not the by-product of a threshold.

## Reproduce

```bash
V=$PWD/.venv/bin/python
bash tools/run_v0416_knee.sh 3                       # ~2 h wall, concurrency 3
$V tools/selection_pool_probe.py --merge <dir>/knee_shard{1,2,3}.json \
    --out-json <dir>/knee_curve_v0416.json
OIN_STRING_EXACT_BOUND=0 $V tools/generator_ab_honest.py --lever OIN_ACCEPT_STRING_EXACT \
    --molecules-file <dir>/pop_wiring40.txt --cohort-dir <cohort> --out-json ab_bound_wiring0.json
```

⚠ Run the A/B tools **from the checkout under test** — `sys.path.insert` beats `PYTHONPATH`.
