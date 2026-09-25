# `measurements/v0.4.19-release-census` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `v0419_rcensus_attach_class_audit.json.gz` | 22660 | `06b9ad0d6efcac83` | tools/census/string_sufficiency.py parseback --side gen / pflags --charge-probe ; tools/attach_class_audit.py   (the per-sweep instruments the attribution reads) |
| `v0419_rcensus_attribution.txt` | 5440 | `8415cb8dfdd3a8ae` | tools/v0419/run_census_release.sh -> tools/census/attribution_table.py --sweep <results-v0.4.19-candidate-sweep> + one path per instrument   (THE CENSUS ON THE v0.4.19 RELEASE SWEEP; C2 = the lane's fix3 canonicality audit, since the encoder changed) |
| `v0419_rcensus_attribution_summary.json` | 4007 | `006e2086d83a74c6` | tools/v0419/run_census_release.sh -> tools/census/attribution_table.py --sweep <results-v0.4.19-candidate-sweep> + one path per instrument   (THE CENSUS ON THE v0.4.19 RELEASE SWEEP; C2 = the lane's fix3 canonicality audit, since the encoder changed) |
| `v0419_rcensus_attribution_table.tsv.gz` | 92681 | `09d7707e4a7436be` | tools/v0419/run_census_release.sh -> tools/census/attribution_table.py --sweep <results-v0.4.19-candidate-sweep> + one path per instrument   (THE CENSUS ON THE v0.4.19 RELEASE SWEEP; C2 = the lane's fix3 canonicality audit, since the encoder changed) |
| `v0419_rcensus_control_on_the_record_sweep.txt` | 5564 | `8ac566db91f21a13` | tools/census/attribution_table.py --write-to <control>   (the parameterised tool must reproduce the frozen v0.4.17 census table byte-for-byte before this release's is believed) |
| `v0419_rcensus_movers_vs_v0418.txt` | 1848 | `722cfc5f4cfdc80b` | tools/census/reattribution_diff.py --old <v0.4.18 release census> --new <this> --movers movers_vs_v0418.txt   (fault transitions against the v0.4.18 release census, split by whether the generated STRUCTURE differs from the v0.4.18 release sweep's) |
| `v0419_rcensus_parseback_gen.jsonl.gz` | 166988 | `9dc1aa74341818fa` | tools/census/string_sufficiency.py parseback --side gen / pflags --charge-probe ; tools/attach_class_audit.py   (the per-sweep instruments the attribution reads) |
| `v0419_rcensus_pflags.jsonl.gz` | 68055 | `b15477135d07096f` | tools/census/string_sufficiency.py parseback --side gen / pflags --charge-probe ; tools/attach_class_audit.py   (the per-sweep instruments the attribution reads) |
| `v0419_rcensus_reattribution_diff.json` | 2411 | `335c650833848d4d` | tools/census/reattribution_diff.py --old <v0.4.18 release census> --new <this> --movers movers_vs_v0418.txt   (fault transitions against the v0.4.18 release census, split by whether the generated STRUCTURE differs from the v0.4.18 release sweep's) |
| `v0419_rcensus_reattribution_diff.txt` | 2922 | `90d7eb2f5277273b` | tools/census/reattribution_diff.py --old <v0.4.18 release census> --new <this> --movers movers_vs_v0418.txt   (fault transitions against the v0.4.18 release census, split by whether the generated STRUCTURE differs from the v0.4.18 release sweep's) |

Source paths at harvest time:

```
v0419_rcensus_attach_class_audit.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_attach_class_audit.json.gz
v0419_rcensus_attribution.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_attribution.txt
v0419_rcensus_attribution_summary.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_attribution_summary.json
v0419_rcensus_attribution_table.tsv.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_attribution_table.tsv.gz
v0419_rcensus_control_on_the_record_sweep.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_control_on_the_record_sweep.txt
v0419_rcensus_movers_vs_v0418.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_movers_vs_v0418.txt
v0419_rcensus_parseback_gen.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_parseback_gen.jsonl.gz
v0419_rcensus_pflags.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_pflags.jsonl.gz
v0419_rcensus_reattribution_diff.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_reattribution_diff.json
v0419_rcensus_reattribution_diff.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.19-release-census/freeze/v0419_rcensus_reattribution_diff.txt
```
