# v0.4.17 baseline sweep -- the sweep OIN_EXACT_DONOR_FOLD's promotion owed

WHY: v0.4.17 changed a shipped DEFAULT (`OIN_EXACT_DONOR_FOLD` joined `levers._DEFAULT_ON`, owner
decision 2026-09-18). That voids the carry-forward licence of the v0.4.14 sweep of record
(77.16% self-consistent / 69.24% verified). This sweep replaces it.

    tools/v0417/launch_sweep.sh  ->  tools/run_sweep.sh <cohort-v0.4.5-5k> <this dir> 6 300
        = tools/test_dataset_roundtrip.py --dataset-dir <cohort> --shard I:6 --mol-timeout 300
          --no-summary        (6 shards, 1-BASED), then tools/rebuild_summary.py
    launched 2026-09-18 12:23 from commit 814abff3 (branch research/v0417-encoder-canonicality),
    as a systemd --user SERVICE unit (OOMPolicy=continue, MemoryMax=14G). 6 h 42 min wall.
    All six shards exited 0. Completeness by REPORT COUNT: 5000 / 5000.

Levers: SHIPPED DEFAULTS. No OIN_* variable was set -- `run_config.json`'s lever block is EMPTY and
a live shard's /proc/<pid>/environ carried none. The exact fold is on only because the registry
says so; that is the thing under test. 12 default-ON levers.
Threads: OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=1, six shards on twelve cores -- like-for-like with
the v0.4.14 sweep of record (its RUN.md explains why: OIN3DGenerator(timeout=) is ADVISORY, so a
contended box shrinks the pool and UNDERSTATES byte_exact).
rdkit 2025.09.3 (pin 2025.9.3). xtb_available=True, optimizer_effective=g-xtb.

## Result  (tools/v0417/post_sweep.sh -> two_numbers.txt; run twice, byte-identical)

    self-consistent   4136 / 5000 = 82.72%      (record 77.16%;  A/B predicted 82.62%)
    VERIFIED          3730 / 5000 = 74.60%      (record 69.24%;  A/B predicted 74.54%)

The reader's control ran first and reproduced the census on the v0.4.14 sweep of record:
5000 / 3858 / 3462.

⚠ `bucket_report.{md,json}` in this directory is the launcher's default SCORED report. The number
is `bucket_report_honest.*`. An inline `smiles_1 == smiles_2_indep` check over-counts (4148).
⚠ `metrics.elapsed_s` is NESTED and SUMS across up to three tiers: a 300 s cap yields ~900 s rows.

## How far to trust it

On the lever's 504 movers this sweep's pass/fail verdict equals the generator A/B's ON arm on
504 / 504. Of the 4496 non-movers 4494 are unchanged; 2 hard_fail timeouts finished this time.
So +278 = +273 (the lever) + 3 (the A/B's own noise-floor flips vs the record) + 2. Book the lever
at +5.46 / +5.30 pts; the other +0.10 / +0.06 is run-to-run, toward fewer timeouts.
smiles_1 equals the exact-fold encoder run on 4988 / 4988. Parse-back re-derived from the NEW
strings moved 1 of 5000 string-fault verdicts.

Full record: docs/agentic-notes/v0.4.17/L1_EXACT_DONOR_FOLD.md
