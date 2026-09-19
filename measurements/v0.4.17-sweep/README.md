# `measurements/v0.4.17-sweep` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `bucket_report_honest.md` | 14566 | `a28060b01ebb25b6` | tools/roundtrip_bucket_report.py --results-dir <dir> |
| `v0417_sweep_RUN.md` | 2723 | `e8ca6fc2a096886f` | tools/v0417/launch_sweep.sh -> tools/run_sweep.sh <cohort-v0.4.5-5k> <out> 6 300   (THE v0.4.17 BASELINE SWEEP: commit 814abff3, SHIPPED DEFAULTS -- the lever block in run_config.json is EMPTY on purpose, that is the thing under test -- 6 shards 1-BASED, --mol-timeout 300, BLAS=1. RUN.md is hand-written provenance) |
| `v0417_sweep_g_verdict.jsonl.gz` | 186462 | `21f1c9c948f00d10` | tools/census/g_vs_input.py --cohort <cohort> --sweep <run> --out <run>   (the census's neutral ruler on the GENERATED STRUCTURES of sweep: the structure half of the VERIFIED predicate. No oinsmiles import) |
| `v0417_sweep_parseback.jsonl.gz` | 156981 | `9bf7cd225108963e` | tools/census/string_sufficiency.py parseback --sweep <sweep> --out <sweep>   (smiles_1 of the NEW sweep read back to a graph with no 3D: the string half of the VERIFIED predicate, RE-DERIVED rather than carried from the census -- 1 of 5000 verdicts moved under the relabeling) |
| `v0417_sweep_run_config.json` | 239 | `cf0c87b57ff2c2f1` | tools/v0417/launch_sweep.sh -> tools/run_sweep.sh <cohort-v0.4.5-5k> <out> 6 300   (THE v0.4.17 BASELINE SWEEP: commit 814abff3, SHIPPED DEFAULTS -- the lever block in run_config.json is EMPTY on purpose, that is the thing under test -- 6 shards 1-BASED, --mol-timeout 300, BLAS=1. RUN.md is hand-written provenance) |
| `v0417_sweep_two_numbers.txt` | 679 | `d0844dd5e35bbce7` | tools/v0417/post_sweep.sh -> tools/v0417/sweep_two_numbers.py --sweep <sweep>   (THE v0.4.17 HEADLINE: self-consistent 4136/5000 = 82.72%, VERIFIED 3730/5000 = 74.60%. The predicate's control runs first and must reproduce the census on the v0.4.14 sweep of record, 5000 / 3858 / 3462, or nothing is printed) |

Source paths at harvest time:

```
bucket_report_honest.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/sweep/bucket_report_honest.md
v0417_sweep_RUN.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/sweep/v0417_sweep_RUN.md
v0417_sweep_g_verdict.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/sweep/v0417_sweep_g_verdict.jsonl.gz
v0417_sweep_parseback.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/sweep/v0417_sweep_parseback.jsonl.gz
v0417_sweep_run_config.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/sweep/v0417_sweep_run_config.json
v0417_sweep_two_numbers.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/sweep/v0417_sweep_two_numbers.txt
```
