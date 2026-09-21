# v0.4.18 baseline sweep -- the sweep the eta construction promotion owed

WHY: v0.4.18 changed two shipped DEFAULTS (`OIN_ETA_COVALENT_TARGET` and `OIN_VDW_EXEMPT_BINDING`
joined `levers._DEFAULT_ON`, owner decision 2026-09-20). That voids the carry-forward licence of
`results-v0.4.17-sweep` (82.72% self-consistent / 74.60% verified). This sweep replaces it.

    tools/v0418/launch_sweep.sh  ->  tools/run_sweep.sh <cohort-v0.4.5-5k> <this dir> 6 300
        = tools/test_dataset_roundtrip.py --dataset-dir <cohort> --shard I:6 --mol-timeout 300
          --no-summary        (6 shards, 1-BASED), then tools/rebuild_summary.py
    launched 2026-09-20 08:48 from commit 2d2708db (branch research/v0418-eta-detached),
    as a systemd --user SERVICE unit (OOMPolicy=continue, MemoryMax=14G). 6 h 24 min wall.
    Completeness by REPORT COUNT: 5000 / 5000.

Levers: SHIPPED DEFAULTS. No OIN_* variable was set -- `run_config.json`'s lever block is EMPTY. The
eta levers are on only because the registry says so; that is the thing under test. 14 default-ON
levers.
Threads: OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=1, six shards on twelve cores -- like-for-like with
the v0.4.17 and v0.4.14 sweeps (OIN3DGenerator(timeout=) is ADVISORY, so a contended box shrinks
the pool and UNDERSTATES byte_exact).

## Result  (tools/v0418/post_sweep.sh -> two_numbers.txt, sweep_vs_ab.txt)

    self-consistent   4305 / 5000 = 86.10%      (v0.4.17 82.72%;  A/B predicted 86.04%)
    VERIFIED          3883 / 5000 = 77.66%      (v0.4.17 74.60%;  A/B predicted 77.60%)

The reader's control ran first and reproduced the census on the v0.4.14 sweep: 5000 / 3858 / 3462.
VERIFIED was computed twice -- with the string faults re-derived from this sweep's own parse-back
(two_numbers.txt) and with them carried from the census table (sweep_vs_ab.txt): 3883 both ways.
⚠ The "(record 77.16% / 69.24%)" labels inside two_numbers.txt are the v0.4.14 record's; the
comparison that matters is the one above.

⚠ `bucket_report.{md,json}` in this directory is the launcher's default SCORED report. The number
is `bucket_report_honest.*`.
⚠ `metrics.elapsed_s` is NESTED and SUMS across up to three tiers: a 300 s cap yields ~900 s rows.

## How far to trust it

The A/B's projection stood on two claims, and sweep_vs_ab.py CHECKED them instead of repeating them:
  * a non-eta molecule cannot reach either lever  -> 3713 / 3713 non-eta structures are
    byte-identical to results-v0.4.17-sweep's (6 more were built in one sweep only)
  * generation is deterministic at this load      -> 1068 / 1068 eta structures are byte-identical
    to the A/B's ON arm's (5 built in one only)
  * smiles_1 differs from the v0.4.17 sweep on 0 of 5000: the levers never touch the encoder.
  * a DEAD promotion would have printed the eta half identical to the OLD sweep; 264 of 1146 are
    (the molecules the levers do not move), not all of them.
Seven verdicts differ from what the projection used, ALL of them rows that were `hard_fail` at the
300 s budget there and finished here (DIFPAM, VUXTOZ, PEDPOG -> byte_exact; PUVWEK, ROMSER, TEQHOM,
WOFGUT -> structural). So +169 self-consistent = +166 (the levers) + 3 (budget boundary), and
+153 verified = +150 + 3. Book the levers at +3.32 / +3.00 pts.
Gains 197 (eta 195) / losses 28 self-consistent; gains 191 (eta 189) / losses 38 VERIFIED -- every
loss is an eta molecule. Runtime: sum elapsed_s 33.18 -> 31.17 h; >30 s 579 -> 515; max 698 -> 660 s.

Full record: docs/agentic-notes/v0.4.18/L2_ETA_DETACHED.md
