# The 5,000-molecule sweep, and what it settles

**Round-trip success 4660/5000 = 93.2 % · `byte_exact` 82.80 %.** Fresh seed-42 cohort from the
25,197 unique tmCAT-tmPHOTO basenames, `--mol-timeout 300`, 6 shards, PASS 1 + PASS 2 complete.
Baseline snapshot: `bucket_report_PASS1_authoritative.{json,md}`.

## ELI5

We generated 3D structures for 5,000 metal complexes and turned each back into a string. 93 % came
back the same. Then we asked *why* the other 7 % didn't — and almost none of it was the notation's
fault: 4 out of 5 failures simply ran out of time before producing anything to compare.

Then we tried the obvious speed-up, and it "recovered" 90 of those failures. It looked like a win.
It wasn't: two-thirds of those 90 were structures whose ligands had fallen off the metal, which the
scoring can't see. The tool we built earlier caught it in the same run.

## 1. Result

| bucket | PASS 1 | final | % |
|---|---|---|---|
| `byte_exact` | 4140 | **4140** | **82.80 %** |
| `key_equal` | 520 | **520** | 10.40 % |
| `facmer_divergent` | 5 | 1 | 0.02 % |
| `structural` | 20 | 9 | 0.18 % |
| `hard_fail` | 300 | 315 | 6.30 % |
| `encode_fail` | 15 | 15 | 0.30 % |

**PASS 2 recovered ZERO molecules.** The g-xTB / FF-reroll recovery ladder moved 15 rows from
"wrong isomer" to "produced nothing" and recovered none. On this cohort it is dead weight — worth
knowing before paying for it again.

## 2. Where the 7 % actually goes — the decomposition that matters

`missed_success_audit`, median failing molecule **300.2 s against a 300 s budget**:

| cause | n | share | tests the notation? |
|---|---|---|---|
| generator timeout | 240 | 70.6 % | **no — never tested** |
| generator produced nothing | 28 | 8.2 % | **no — never tested** |
| canonicalization noise | 32 | 9.4 % | test too strict |
| output names a different isomer | 25 | 7.4 % | ambiguous |
| encoder refused the input | 15 | 4.4 % | encoder coverage |

**78.8 % of failures never test the notation.** The notation-attributable gap is **57/5000 ≈ 1.1 %** —
where the notation is actually exercised it is ~98.9 % correct.

**Consequence: the `<30 s` goal and the `100 %` goal are ONE goal.** Both are bounded by generator
compute, not by the notation. It also means every pass-rate figure in this project is partly a
runtime measurement, which is why v0.4.4's 11 "regressions" were all timeouts.

## 3. The throughput path, tested — and why it must not ship default-ON

If failures are timeouts, a faster acceptance predicate should convert them into measurements. Tested
directly: the 340 failing molecules re-run with `OIN_ACCEPT_SCORED=1` at a **60 s** cap — deliberately
5× harsher than the sweep's 300 s, so every pass is a conservative lower bound.

| | |
|---|---|
| apparent recoveries | **90 / 340 = 26.5 %** |
| median elapsed of those | 13.7 s (74 of 90 **under 30 s**) |
| apparent corpus rate | 4660 + 90 = **95.0 %** |

That is where the story would have ended without an instrument. `report["coordination"]` is now
always-on, so the same run reports what those structures actually are:

| coordination verdict of the 90 | n |
|---|---|
| **DEGRADED — ligands came off the metal** | **60 (66.7 %)** |
| BOUNDARY — held within 0.1 Å of the cutoff | 21 |
| **INTACT — genuine recovery** | **9** |

**So +90 → 95.0 % is really +9 → 93.4 %.** Two-thirds of the accuracy the lever appears to buy is
detached ligands scored as passes. This confirms `ACCEPT_SCORED_v0.4.7.md`'s G2 gate first-hand and
at scale, and it is the decisive argument against promoting the lever globally: **flipping it would
manufacture 60 phantom passes that no existing sweep or gate could detect.**

The runtime win is real and remains worth having under scope — throughput and generator-benchmarking
runs, where the pass column is not the deliverable. Data: `docs/ACCEPT_SCORED_RECOVERY_ARM.json`.

## 4. What this says about reaching 100 %

1. **The notation is ~98.9 % correct where tested.** That clears a 97-98 % bar. Remaining notation
   work is 57 molecules, itemised in the audit's informative remainder.
2. **The residual 6.3 % is generator throughput and capability**, and it is not reachable by encoder
   or notation work. Boron is the clearest case: 0/10 cage molecules assemble at all.
3. **Speed alone does not buy accuracy honestly.** The one available 3–4× speed-up converts
   timeouts into *wrong structures* two times out of three. Any future throughput work must be
   scored with `coordination` (and ideally `OIN_INDEP_SCORE`) or it will report progress it has not
   made.
4. **On an honest metric the headline goes DOWN, not up** — the false-positive audit puts the
   over-statement at ~5.7 points net. 93.2 % scored is roughly 88 % honest. The instrument to
   measure that now exists and is off by default only because it costs a second encode.

## 5. Reproduce

```bash
export PYTHONPATH=$PWD/src; V=.venv/bin/python
D=$PWD/tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep
$V tools/roundtrip_bucket_report.py --results-dir "$D"
$V -m tools.injectivity.missed_success_audit --sweep "$D"    # TIMEOUT_S=300 must match the budget

# the throughput arm, on the failing population only
OIN_ACCEPT_SCORED=1 $V tools/test_dataset_roundtrip.py --dataset-dir <340-failures> \
    --output-dir results-v047-acceptscored-recovery --shard 1:4 --no-summary --mol-timeout 60
```
