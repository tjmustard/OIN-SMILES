# OIN-SMILES v0.4.4 — round-trip accuracy (regression sweep + full-set projection)

What the v0.4.4 default changed for XYZ ↔ OIN round-trip fidelity, measured on a clean
single-commit regression sweep, plus an honest projection to the full dataset. Companion to the
`CHANGELOG.md` `[0.4.4]` entry (mechanisms), `docs/DIRECT_DG_VALIDATION.md` (the direct-assembly
decision), and `docs/KNOWN_LIMITATIONS.md` (what is out of scope). Raw per-molecule data:
`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.4-regression/REGRESSION_REPORT.md` (results dir gitignored).

## Baseline

The v0.4.0 continuous accumulator processed **25,197** molecules of the ~26,230-file tmCAT/tmPHOTO
set: **22,280 pass / 2,917 fail = 88.4%** (quick mode, 30 s, older comparison). This is the
reference the v0.4.4 sweep measures against.

## What was measured — the regression sweep

A full-quality sweep of the v0.4.4 default over a **3,917-molecule** set: **all 2,917 v0.4.0
failures** + a **1,000-molecule seed-42 sample of the 22,280 v0.4.0 successes**. Single commit
`f50a199e` (the v0.4.4 release tip; the stack has since been rebased with identical content), full
quality (`timeout=300`s, real g-xTB, `seed=42`). Both the v0.4.0 and v0.4.4 results were scored
with the **same** fac/mer `oin.compare` key, so key-version and rdkit-version drift cannot
fabricate a regression. `OK = byte_exact | key_equal`.

**Headline: 11 regressions, 1,092 fixes, net +1,081 round-trip-OK.**

| cohort | n | v0.4.0 OK | v0.4.4 OK | Δ | basis |
|---|---:|---|---|---:|---|
| previously passing (guard) | 1000 (of 22,280) | 999 (99.9%) | 988 (98.8%) | −11 | random sample |
| previously failing (target) | 2917 | 0 | 1092 (37.4%) | +1092 | full cohort, measured |
| sweep total | 3917 | 999 (25.5%) | 2080 (53.1%) | +1081 | — |

**All 11 regressions are 300 s generation timeouts** (10 UFF, 1 g-xTB), **zero** correctness or
notation regressions. Each is a 69–157-atom molecule that round-tripped in 6–29 s under v0.4.0
quick mode but exceeded the 300 s cap under full quality (full pool + direct-DG). This is a
compute-time artifact of the config asymmetry, not a wrong answer — see
`docs/KNOWN_LIMITATIONS.md`.

### Bucket movement (on the 3,917 set)

| bucket | v0.4.0 | v0.4.4 |
|---|---:|---:|
| byte_exact | 851 | 1759 |
| key_equal | 148 | 321 |
| facmer_divergent | 36 | 2 |
| structural | 349 | 42 |
| hard_fail | 2116 | 1552 |
| encode_fail | 417 | 241 |

The `structural` bucket collapsed 349 → 42 and `encode_fail` 417 → 241; the fixes come from
`hard_fail`/`structural`/`encode_fail`/`facmer_divergent` → `byte_exact`/`key_equal`.

## Full-set pass rate — a projection, not a measurement

The full 25,197-molecule set was **not** re-run under v0.4.4. The number below extrapolates the two
measured cohorts back onto the v0.4.0 universe.

| cohort | size | v0.4.4 pass rate | projected passes |
|---|---:|---|---:|
| previously passing | 22,280 | 98.8% (sample) | ~22,010 |
| previously failing | 2,917 | 37.4% (measured) | 1,092 |
| **total** | **25,197** | **~91.7%** | **~23,100** |

- **Point estimate: ~91.7%**, sampling-noise range **91.1–92.3%** (the success piece is a
  4.5% sample; the failure piece is exact).
- **Correctness upper bound ~92.8%** — the 11/1000 success-cohort regressions are all 300 s
  timeouts that would pass at quick config or with more compute; exclude them and the projection is
  (22,280 + 1,092) / 25,197.
- **Net vs v0.4.0's 88.4%: roughly +3 to +4 points.**

State it as **~92%, ±~1, projected** — not a measured headline.

## Why a single headline number is fraught (as in v0.4.2)

1. **The measurement lens changed.** This is the fac/mer-key round-trip rate under full quality
   (300 s); the 88.4% baseline was quick mode (30 s) with the older comparison. The ~+3-point gain
   mixes a real accuracy improvement with a stricter, cleaner metric.
2. **The success-cohort drag is a compute artifact, not lost correctness.** The extrapolated ~245
   timeout regressions are slow big molecules, not wrong answers.
3. **It is a projection.** The failure cohort was measured in full; the success cohort was sampled
   at ~4.5%.

## To replace the estimate with a measurement

Run the full 25,197 set under the v0.4.4 default, same harness as the regression sweep
(`tools/test_dataset_roundtrip.py` → `rebuild_summary.py` → `tools/roundtrip_bucket_report.py`).
**Quick mode** (`--quick`, 30 s) is apples-to-apples with the 88.4% baseline and far faster than
the multi-day full-quality tail; **full quality** (300 s) matches the sweep above but re-incurs the
timeout tail. Prefer quick mode for a methodology-matched headline.

## Reproducing the sweep

1. Manifest: all `status=="failed"` (2,917) + `random.seed(42); random.sample(sorted(successes),
   1000)` from `results-v0.4.0/summary_roundtrip.json` → 3,917 names.
2. Deduped symlink inputs (`regression_inputs/`, exact v0.4.0 copies), 6-shard `--no-summary`
   sweep at `--mol-timeout 300` → `tools/rebuild_summary.py`.
3. `tools/roundtrip_bucket_report.py` on the v0.4.4 results **and** on `results-v0.4.0` restricted
   with `--only` to the same 3,917; join per molecule by bucket transition.
