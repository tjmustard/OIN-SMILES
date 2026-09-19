# `measurements/v0.4.18-census` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `v0418_census_attach_class_audit.json.gz` | 22823 | `5e400050cd348c29` | tools/census/string_sufficiency.py parseback --side gen | pflags --charge-probe ; tools/attach_class_audit.py   (the three SWEEP-DEPENDENT instruments, derived on the v0.4.17 sweep. The other inputs are frozen where they were made: g_verdict + parseback in v0.4.17-sweep/, e_selfconsistency_exact in v0.4.17/, the mirror ruler + twin clusters in v0.4.17-census/) |
| `v0418_census_attribution.txt` | 5501 | `511b722d1714d5e4` | tools/census/attribution_table.py --sweep <results-v0.4.17-sweep> --out <reattr> + one path per instrument   (THE CENSUS RE-RUN ON THE v0.4.17 SWEEP: 5000 rows, one fault each, UNATTRIBUTED 0. FAIL 864 = 17.28 pts: G 9.62 / E2 3.50 / E1 3.38 / P 0.58 / data 0.20. NONE = 3730 = the VERIFIED count sweep_two_numbers.py printed, or the tool aborts) |
| `v0418_census_attribution_summary.json` | 4248 | `889e03a12c28e1f6` | tools/census/attribution_table.py --sweep <results-v0.4.17-sweep> --out <reattr> + one path per instrument   (THE CENSUS RE-RUN ON THE v0.4.17 SWEEP: 5000 rows, one fault each, UNATTRIBUTED 0. FAIL 864 = 17.28 pts: G 9.62 / E2 3.50 / E1 3.38 / P 0.58 / data 0.20. NONE = 3730 = the VERIFIED count sweep_two_numbers.py printed, or the tool aborts) |
| `v0418_census_attribution_table.tsv.gz` | 94386 | `06a31ad75da36118` | tools/census/attribution_table.py --sweep <results-v0.4.17-sweep> --out <reattr> + one path per instrument   (THE CENSUS RE-RUN ON THE v0.4.17 SWEEP: 5000 rows, one fault each, UNATTRIBUTED 0. FAIL 864 = 17.28 pts: G 9.62 / E2 3.50 / E1 3.38 / P 0.58 / data 0.20. NONE = 3730 = the VERIFIED count sweep_two_numbers.py printed, or the tool aborts) |
| `v0418_census_control_on_the_record_sweep.txt` | 5563 | `99e4b5b42dd0e388` | tools/census/attribution_table.py --write-to <scratch>   (CONTROL, run first: the parameterised tool on the v0.4.14 record sweep reproduces measurements/v0.4.17-census/attribution_table.tsv.gz BYTE-FOR-BYTE) |
| `v0418_census_parseback_gen.jsonl.gz` | 166656 | `6d003af783390484` | tools/census/string_sufficiency.py parseback --side gen | pflags --charge-probe ; tools/attach_class_audit.py   (the three SWEEP-DEPENDENT instruments, derived on the v0.4.17 sweep. The other inputs are frozen where they were made: g_verdict + parseback in v0.4.17-sweep/, e_selfconsistency_exact in v0.4.17/, the mirror ruler + twin clusters in v0.4.17-census/) |
| `v0418_census_pflags.jsonl.gz` | 68159 | `d68811487b7a3a84` | tools/census/string_sufficiency.py parseback --side gen | pflags --charge-probe ; tools/attach_class_audit.py   (the three SWEEP-DEPENDENT instruments, derived on the v0.4.17 sweep. The other inputs are frozen where they were made: g_verdict + parseback in v0.4.17-sweep/, e_selfconsistency_exact in v0.4.17/, the mirror ruler + twin clusters in v0.4.17-census/) |
| `v0418_census_reattribution_diff.json` | 2335 | `3007942b6b0b1f4b` | tools/census/reattribution_diff.py --old <census table> --new <new table> --movers <504>   (where each fault WENT, split by mover: 4482 of 4496 non-movers keep theirs -- the noise floor of attribution is 14 molecules / 0.28 pts) |
| `v0418_census_reattribution_diff.txt` | 2973 | `fa9fbfa59346857a` | tools/census/reattribution_diff.py --old <census table> --new <new table> --movers <504>   (where each fault WENT, split by mover: 4482 of 4496 non-movers keep theirs -- the noise floor of attribution is 14 molecules / 0.28 pts) |

Source paths at harvest time:

```
v0418_census_attach_class_audit.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_attach_class_audit.json.gz
v0418_census_attribution.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_attribution.txt
v0418_census_attribution_summary.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_attribution_summary.json
v0418_census_attribution_table.tsv.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_attribution_table.tsv.gz
v0418_census_control_on_the_record_sweep.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_control_on_the_record_sweep.txt
v0418_census_parseback_gen.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_parseback_gen.jsonl.gz
v0418_census_pflags.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_pflags.jsonl.gz
v0418_census_reattribution_diff.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_reattribution_diff.json
v0418_census_reattribution_diff.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.17-reattribution/freeze/v0418_census_reattribution_diff.txt
```
