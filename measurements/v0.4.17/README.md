# `measurements/v0.4.17` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `v0417_ab_gains.txt` | 3864 | `ee391413d04a1ef5` | tools/v0417/run_generator_ab.sh <union movers> ; tools/v0417/post_generator_ab.sh   (L1b: the REAL sweep harness over the 504 movers, OFF and ON arms run SIMULTANEOUSLY, commit de4a9d2c. Self-consistent 276 gains / 3 losses, VERIFIED 269 / 4; dead-lever check 392; noise floor 3/504 vs the sweep of record) |
| `v0417_ab_losses.txt` | 42 | `a2002ec4941735e2` | tools/v0417/run_generator_ab.sh <union movers> ; tools/v0417/post_generator_ab.sh   (L1b: the REAL sweep harness over the 504 movers, OFF and ON arms run SIMULTANEOUSLY, commit de4a9d2c. Self-consistent 276 gains / 3 losses, VERIFIED 269 / 4; dead-lever check 392; noise floor 3/504 vs the sweep of record) |
| `v0417_ab_off_bucket_report_honest.json.gz` | 27532 | `2a44e11c54839138` | tools/v0417/run_generator_ab.sh <union movers> ; tools/v0417/post_generator_ab.sh   (L1b: the REAL sweep harness over the 504 movers, OFF and ON arms run SIMULTANEOUSLY, commit de4a9d2c. Self-consistent 276 gains / 3 losses, VERIFIED 269 / 4; dead-lever check 392; noise floor 3/504 vs the sweep of record) |
| `v0417_ab_off_g_verdict.jsonl.gz` | 19129 | `d19a15fff7334528` | tools/census/g_vs_input.py --cohort <cohort> --sweep <run> --out <run>   (the census's neutral ruler on the GENERATED STRUCTURES of ab_off: the structure half of the VERIFIED predicate. No oinsmiles import) |
| `v0417_ab_on_bucket_report_honest.json.gz` | 25695 | `2a8322c0252e1a0d` | tools/v0417/run_generator_ab.sh <union movers> ; tools/v0417/post_generator_ab.sh   (L1b: the REAL sweep harness over the 504 movers, OFF and ON arms run SIMULTANEOUSLY, commit de4a9d2c. Self-consistent 276 gains / 3 losses, VERIFIED 269 / 4; dead-lever check 392; noise floor 3/504 vs the sweep of record) |
| `v0417_ab_on_g_verdict.jsonl.gz` | 20023 | `11725c145bc9b276` | tools/census/g_vs_input.py --cohort <cohort> --sweep <run> --out <run>   (the census's neutral ruler on the GENERATED STRUCTURES of ab_on: the structure half of the VERIFIED predicate. No oinsmiles import) |
| `v0417_autofold_audit.json.gz` | 114066 | `149ecb219e714b42` | tools/v0417/autofold_audit.py   (L1a OFFLINE GATE, exact because the slot post-pass is a pure string function: bucket fold vs exact fold applied to the census's stored rotation-only strings. Positive control 5792/5792; 550 of the 568 vetoed achiral pairs unified; 18 kept split are chiral by wrap; 17 of the 36 protected pairs unified) |
| `v0417_automorphism_extension_check.json` | 277 | `af9aa60a7b7d930d` | tools/v0417/automorphism_extension_check.py   (does pendant stripping invent a symmetry? 13716 pruned automorphisms over 7968 fragment checks, 0 without a full-graph extension; negative control refused 7953/7953) |
| `v0417_e_selfconsistency_exact.jsonl.gz` | 252493 | `9dce5618ff1f4889` | OIN_EXACT_DONOR_FOLD=1 tools/census/e_selfconsistency.py --name e_selfconsistency_exact   (census C2 re-run unchanged under the lever, same seeds: eleven presentations of every input. Diff against census_e_selfconsistency.jsonl.gz: achiral-and-differs 670 -> 120, slot-level renumber drift 494 -> 33, KEY-level 255 -> 255 with membership changed on 0) |
| `v0417_generated_side_movers.jsonl.gz` | 193968 | `b401ae5817d80e44` | tools/v0417/generated_side_movers.py   (every stored generated structure of the v0.4.14 sweep encoded lever OFF and ON: 286 of 4743 move. OFF == the sweep of record on 4743/4743) |
| `v0417_generator_ab_report.txt` | 2253 | `2177d64c7fbf828d` | tools/v0417/run_generator_ab.sh <union movers> ; tools/v0417/post_generator_ab.sh   (L1b: the REAL sweep harness over the 504 movers, OFF and ON arms run SIMULTANEOUSLY, commit de4a9d2c. Self-consistent 276 gains / 3 losses, VERIFIED 269 / 4; dead-lever check 392; noise floor 3/504 vs the sweep of record) |
| `v0417_pop_ab_union_movers.txt` | 7056 | `6cfb4bcd5387a4a9` | tools/v0417/exact_fold_live_report.py (input side, 392) / union with the generated side   (SAMPLE MEMBERSHIP of the generator A/B: 504 = 392 u 286, 174 in both. Every other molecule gets the same two strings in both arms and is unchanged by construction -- which the sweep then measured: 4494 of 4496) |
| `v0417_pop_generated_side_movers.txt` | 4004 | `389a3d4d7e17c797` | tools/v0417/generated_side_movers.py   (every stored generated structure of the v0.4.14 sweep encoded lever OFF and ON: 286 of 4743 move. OFF == the sweep of record on 4743/4743) |
| `v0417_pop_input_side_movers.txt` | 5488 | `24475d6d3280712a` | tools/v0417/exact_fold_live_report.py (input side, 392) / union with the generated side   (SAMPLE MEMBERSHIP of the generator A/B: 504 = 392 u 286, 174 in both. Every other molecule gets the same two strings in both arms and is unchanged by construction -- which the sweep then measured: 4494 of 4496) |

Source paths at harvest time:

```
v0417_ab_gains.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_ab_gains.txt
v0417_ab_losses.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_ab_losses.txt
v0417_ab_off_bucket_report_honest.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_ab_off_bucket_report_honest.json.gz
v0417_ab_off_g_verdict.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_ab_off_g_verdict.jsonl.gz
v0417_ab_on_bucket_report_honest.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_ab_on_bucket_report_honest.json.gz
v0417_ab_on_g_verdict.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_ab_on_g_verdict.jsonl.gz
v0417_autofold_audit.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_autofold_audit.json.gz
v0417_automorphism_extension_check.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_automorphism_extension_check.json
v0417_e_selfconsistency_exact.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_e_selfconsistency_exact.jsonl.gz
v0417_generated_side_movers.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_generated_side_movers.jsonl.gz
v0417_generator_ab_report.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_generator_ab_report.txt
v0417_pop_ab_union_movers.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_pop_ab_union_movers.txt
v0417_pop_generated_side_movers.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_pop_generated_side_movers.txt
v0417_pop_input_side_movers.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-exactfold/freeze/l1/v0417_pop_input_side_movers.txt
```
