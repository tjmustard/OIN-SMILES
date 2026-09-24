# v0.4.18 RELEASE sweep -- what ships: both eta promotions as shipped defaults

WHY: v0.4.18 changed three shipped DEFAULTS in two steps -- `OIN_ETA_COVALENT_TARGET` +
`OIN_VDW_EXEMPT_BINDING` (owner decision 2026-09-20; measured by results-v0.4.18-sweep, 86.10 /
77.66) and `OIN_ETA_RETARGET` (owner delegation 2026-09-23). The second voids the first sweep's
carry-forward licence. This sweep is the baseline of record for v0.4.18.

    SWEEP_TAG=v0.4.18-release tools/v0418/launch_sweep.sh
        -> tools/run_sweep.sh <cohort-v0.4.5-5k> <this dir> 6 300
        = tools/test_dataset_roundtrip.py --dataset-dir <cohort> --shard I:6 --mol-timeout 300
          --no-summary        (6 shards, 1-BASED), then tools/rebuild_summary.py
    launched 2026-09-23 10:52 from commit cdb76eca (branch research/v0418-eta-selection),
    as a systemd --user SERVICE unit (OOMPolicy=continue, MemoryMax=14G). 6 h 29 min wall.
    Completeness by REPORT COUNT: 5000 / 5000.

Levers: SHIPPED DEFAULTS. No OIN_* variable was set -- `run_config.json`'s lever block is EMPTY.
15 default-ON levers. Threads: OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS=1, six shards on twelve cores
-- like-for-like with every sweep of record since v0.4.14.

## Result  (SWEEP_TAG=v0.4.18-release tools/v0418/post_sweep.sh -> two_numbers.txt, sweep_vs_ab.txt)

    self-consistent   4347 / 5000 = 86.94%      (v0.4.18-sweep 86.10%;  A/B predicted 87.00%)
    VERIFIED          3920 / 5000 = 78.40%      (v0.4.18-sweep 77.66%;  A/B predicted 78.42%)
    v0.4.17 sweep of record, for the release delta: 82.72% / 74.60%  ->  +4.22 / +3.80 pts

The reader's control ran first and reproduced the census on the v0.4.14 sweep: 5000 / 3858 / 3462.
VERIFIED computed two ways (string faults re-derived from this sweep's parse-back; carried from the
census table): 3920 both. ⚠ two_numbers.txt's "(record ...)" labels are the v0.4.14 record's.
⚠ `bucket_report.{md,json}` is the launcher's default SCORED report; the number is
`bucket_report_honest.*`. ⚠ `metrics.elapsed_s` is NESTED and a SUM over up to three tiers.

## How far to trust it

sweep_vs_ab.py CHECKED the two claims the selection A/B's projection stood on:
  * a non-eta molecule cannot reach the lever   -> 3713 / 3713 non-eta structures byte-identical to
    results-v0.4.18-sweep (and so to results-v0.4.17-sweep); 7 built in one sweep only
  * generation is deterministic at this load    -> 1070 / 1070 eta structures byte-identical to the
    A/B's retarget arm (results-v0.4.18-eta-selection/ab_retarget); 8 built in one only
  * smiles_1 differs from the v0.4.17 sweep on 0 of 5000: no promotion in v0.4.18 touches the encoder
  * 912 of 1146 eta structures are identical to results-v0.4.18-sweep: the molecules OIN_ETA_RETARGET
    does not move (the A/B said 158 of 1067 differ).
Eleven verdicts differ from what the projection used and ALL are rows at the 300 s budget
(hard_fail <-> structural / byte_exact, elapsed ~300 s): ECIGAZ, NOEPOR, PEDPEW, QEBKUG, QOCHAT,
ROMSER, TEQHOM, VUXTOZ, WOFGUT, XOSCIT, YOQMAT. So 4347 / 3920 = the A/B's 4350 / 3921 - 3 / - 1 of
timing. Gains 47 (eta 47) / losses 5 self-consistent; 43 (eta 43) / 6 verified, against
results-v0.4.18-sweep. Buckets: structural 320 -> 284, key_equal 129 -> 119, hard_fail 228 -> 230.
Runtime: sum elapsed_s 31.17 -> 31.95 h; >30 s 515 -> 521; max 660 -> 628 s.

Full record: docs/agentic-notes/v0.4.18/SELECTION_ETA_RETARGET.md (section 6) and L2_ETA_DETACHED.md
