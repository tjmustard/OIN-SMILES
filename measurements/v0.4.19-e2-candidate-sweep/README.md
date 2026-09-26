# `measurements/v0.4.19-e2-candidate-sweep` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `v0419_e2sweep_bucket_report_honest.json.gz` | 224267 | `88c3e066a1b0baf7` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS=... SWEEP_SHIPPED=... SWEEP_PROJECTION=e2f tools/v0419/post_candidate_sweep.sh   (honest buckets, the neutral ruler, parse-back under the levers and under the shipped reader; tools/v0417/sweep_two_numbers.py) |
| `v0419_e2sweep_bucket_report_honest.md` | 8671 | `38312b21f7747a74` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS=... SWEEP_SHIPPED=... SWEEP_PROJECTION=e2f tools/v0419/post_candidate_sweep.sh   (honest buckets, the neutral ruler, parse-back under the levers and under the shipped reader; tools/v0417/sweep_two_numbers.py) |
| `v0419_e2sweep_g_verdict.jsonl.gz` | 190509 | `ec0e3b2b866e311a` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS=... SWEEP_SHIPPED=... SWEEP_PROJECTION=e2f tools/v0419/post_candidate_sweep.sh   (honest buckets, the neutral ruler, parse-back under the levers and under the shipped reader; tools/v0417/sweep_two_numbers.py) |
| `v0419_e2sweep_parseback.jsonl.gz` | 155682 | `60aa61ec2b9f2978` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS=... SWEEP_SHIPPED=... SWEEP_PROJECTION=e2f tools/v0419/post_candidate_sweep.sh   (honest buckets, the neutral ruler, parse-back under the levers and under the shipped reader; tools/v0417/sweep_two_numbers.py) |
| `v0419_e2sweep_parseback_gen.jsonl.gz` | 166806 | `7e7fb9c120b7aa02` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS=... SWEEP_SHIPPED=... SWEEP_PROJECTION=e2f tools/v0419/post_candidate_sweep.sh   (honest buckets, the neutral ruler, parse-back under the levers and under the shipped reader; tools/v0417/sweep_two_numbers.py) |
| `v0419_e2sweep_parseback_shippedreader.jsonl.gz` | 155682 | `60aa61ec2b9f2978` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS=... SWEEP_SHIPPED=... SWEEP_PROJECTION=e2f tools/v0419/post_candidate_sweep.sh   (honest buckets, the neutral ruler, parse-back under the levers and under the shipped reader; tools/v0417/sweep_two_numbers.py) |
| `v0419_e2sweep_run_config.json` | 302 | `7751ff74a739c6bb` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS='-E OIN_N_VALENCE_2=1 -E OIN_CANONICAL_RESONANCE=1' tools/v0419/launch_candidate_sweep.sh   (the full 5,000, 6 shards, --mol-timeout 300; the levers block must be exactly the two) |
| `v0419_e2sweep_two_numbers.txt` | 691 | `53ff36bea5b246e3` | SWEEP_TAG=v0.4.19-e2-candidate SWEEP_LEVERS=... SWEEP_SHIPPED=... SWEEP_PROJECTION=e2f tools/v0419/post_candidate_sweep.sh   (honest buckets, the neutral ruler, parse-back under the levers and under the shipped reader; tools/v0417/sweep_two_numbers.py) |
| `v0419_e2sweep_unchanged_string_movers.json` | 7485 | `6120eed788ae06b9` | tools/v0419/candidate_vs_projection.py --projection e2f   (every row where the sweep and the e2f projection disagree, with elapsed_s in both runs; the generator-reach table) |
| `v0419_e2sweep_vs_projection.json` | 5397 | `a86574044f10afd3` | tools/v0419/candidate_vs_projection.py --projection e2f   (every row where the sweep and the e2f projection disagree, with elapsed_s in both runs; the generator-reach table) |
| `v0419_e2sweep_vs_projection.txt` | 4701 | `13ebd46196b25954` | tools/v0419/candidate_vs_projection.py --projection e2f   (every row where the sweep and the e2f projection disagree, with elapsed_s in both runs; the generator-reach table) |

Source paths at harvest time:

```
v0419_e2sweep_bucket_report_honest.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_bucket_report_honest.json.gz
v0419_e2sweep_bucket_report_honest.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_bucket_report_honest.md
v0419_e2sweep_g_verdict.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_g_verdict.jsonl.gz
v0419_e2sweep_parseback.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_parseback.jsonl.gz
v0419_e2sweep_parseback_gen.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_parseback_gen.jsonl.gz
v0419_e2sweep_parseback_shippedreader.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_parseback_shippedreader.jsonl.gz
v0419_e2sweep_run_config.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_run_config.json
v0419_e2sweep_two_numbers.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_two_numbers.txt
v0419_e2sweep_unchanged_string_movers.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_unchanged_string_movers.json
v0419_e2sweep_vs_projection.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_vs_projection.json
v0419_e2sweep_vs_projection.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-e2-candidate-sweep/freeze/v0419_e2sweep_vs_projection.txt
```
