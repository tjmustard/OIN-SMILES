# `measurements/v0.4.20-release-census` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `v0420_rcensus_attach_class_audit.json.gz` | 22632 | `5a5e4c84a68f4ea4` | tools/census/string_sufficiency.py parseback --side gen / pflags --charge-probe ; tools/attach_class_audit.py   (the per-sweep instruments the attribution reads) |
| `v0420_rcensus_attribution.txt` | 5436 | `c49aa755f6390a2c` | tools/v0420/run_census_release.sh -> tools/census/attribution_table.py --sweep <results-v0.4.19-e2-candidate-sweep> + one path per instrument   (THE CENSUS ON THE v0.4.20 RELEASE SWEEP; C2 = the E2 lane's e2f canonicality audit, the same two levers) |
| `v0420_rcensus_attribution_summary.json` | 4002 | `31ebe9d21bd9d221` | tools/v0420/run_census_release.sh -> tools/census/attribution_table.py --sweep <results-v0.4.19-e2-candidate-sweep> + one path per instrument   (THE CENSUS ON THE v0.4.20 RELEASE SWEEP; C2 = the E2 lane's e2f canonicality audit, the same two levers) |
| `v0420_rcensus_attribution_table.tsv.gz` | 91917 | `bf997bc14fef667a` | tools/v0420/run_census_release.sh -> tools/census/attribution_table.py --sweep <results-v0.4.19-e2-candidate-sweep> + one path per instrument   (THE CENSUS ON THE v0.4.20 RELEASE SWEEP; C2 = the E2 lane's e2f canonicality audit, the same two levers) |
| `v0420_rcensus_control_on_the_record_sweep.txt` | 5564 | `ab4a0d7552a20d14` | tools/census/attribution_table.py --write-to <control>   (the parameterised tool must reproduce the frozen v0.4.17 census table byte-for-byte before this release's is believed) |
| `v0420_rcensus_movers_vs_v0419.txt` | 2016 | `46e4c76aa628cf3d` | tools/census/reattribution_diff.py --old <v0.4.19 release census> --new <this> --movers movers_vs_v0419.txt   (fault transitions against the v0.4.19 release census, split by whether the generated STRUCTURE differs from the v0.4.19 release sweep's) |
| `v0420_rcensus_parseback_gen.jsonl.gz` | 166806 | `7e7fb9c120b7aa02` | tools/census/string_sufficiency.py parseback --side gen / pflags --charge-probe ; tools/attach_class_audit.py   (the per-sweep instruments the attribution reads) |
| `v0420_rcensus_pflags.jsonl.gz` | 67896 | `96b60b25cb0fa27b` | tools/census/string_sufficiency.py parseback --side gen / pflags --charge-probe ; tools/attach_class_audit.py   (the per-sweep instruments the attribution reads) |
| `v0420_rcensus_reattribution_diff.json` | 2086 | `4835348d4c40adfd` | tools/census/reattribution_diff.py --old <v0.4.19 release census> --new <this> --movers movers_vs_v0419.txt   (fault transitions against the v0.4.19 release census, split by whether the generated STRUCTURE differs from the v0.4.19 release sweep's) |
| `v0420_rcensus_reattribution_diff.txt` | 2818 | `d187ffa890199c26` | tools/census/reattribution_diff.py --old <v0.4.19 release census> --new <this> --movers movers_vs_v0419.txt   (fault transitions against the v0.4.19 release census, split by whether the generated STRUCTURE differs from the v0.4.19 release sweep's) |

Source paths at harvest time:

```
v0420_rcensus_attach_class_audit.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_attach_class_audit.json.gz
v0420_rcensus_attribution.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_attribution.txt
v0420_rcensus_attribution_summary.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_attribution_summary.json
v0420_rcensus_attribution_table.tsv.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_attribution_table.tsv.gz
v0420_rcensus_control_on_the_record_sweep.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_control_on_the_record_sweep.txt
v0420_rcensus_movers_vs_v0419.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_movers_vs_v0419.txt
v0420_rcensus_parseback_gen.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_parseback_gen.jsonl.gz
v0420_rcensus_pflags.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_pflags.jsonl.gz
v0420_rcensus_reattribution_diff.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_reattribution_diff.json
v0420_rcensus_reattribution_diff.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.20-release-census/freeze/v0420_rcensus_reattribution_diff.txt
```
