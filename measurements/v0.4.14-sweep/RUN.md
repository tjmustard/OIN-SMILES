# v0.4.14 baseline sweep

WHY: every byte_exact figure since v0.4.6 is an offline re-score of v0.4.6-era generated
structures. v0.4.9-v0.4.14 all changed what the generator builds, so the ABSOLUTE has
drifted even though each release's DELTA was measured. This sweep re-establishes it.

    tools/test_dataset_roundtrip.py --dataset-dir <cohort-v0.4.5-5k> --shard I:6
        --mol-timeout 300 --no-summary        (6 shards, 1-BASED)
    then tools/rebuild_summary.py

Levers: SHIPPED DEFAULTS, all explicitly unset before launch (11 default-ON incl.
OIN_RESONANCE_DONOR_FOLD and OIN_INDEP_SCORE).
⚠ --shard is 1-BASED. `--shard 0:6` exits 2, and launching 0..5 silently drops 1/6.
⚠ metrics.elapsed_s is NESTED and SUMS across up to three tiers: a 300s cap yields ~900s rows.

## Launch note (2026-07-29)
First attempt used `systemd-run --user --scope -p OOMPolicy=continue -p MemoryMax=14G` and
failed: **`OOMPolicy` is a SERVICE property and `--scope` rejects it** ("Unknown assignment").
The recorded precedent for that flag was a service unit, not a scope.

Dropped systemd entirely. The six shards are independent processes, so an OOM kill takes one
shard rather than the run -- which is what `OOMPolicy=continue` was buying. Completeness is
verified by REPORT COUNT (must reach 5000), not by exit status, so a silently short shard
cannot pass as a finished sweep.

## Thread caps (2026-07-29, restart #2)
First run at 6 shards drove load to ~30 on 12 cores: each molecule's child was taking
125-235% CPU from BLAS, so six shards oversubscribed badly.

That is not merely slow, and the reason it forced a restart is specific:
**`OIN3DGenerator(timeout=)` is ADVISORY** -- the embed loop checks its deadline BETWEEN
attempts. Under CPU starvation fewer attempts fit inside the budget, the pool is smaller, and
the sweep would **understate** byte_exact. A contended baseline is a biased baseline.

So: OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=1, six shards ~= six cores on twelve.
⚠ This DEVIATES from however the v0.4.6 sweep was run, so per-molecule elapsed_s is not
strictly like-for-like with the 994/5000 figure. Accuracy is the point of this run; the
runtime comparison inherits a caveat and is stated with one.
