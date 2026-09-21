# `measurements/v0.4.18-selection` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `v0418_sel_ab_commits.tsv` | 289 | `96bfe5b8c671e6bd` | tools/v0418/freeze_stage_selection.py   (MADE at staging time: per molecule and arm the sha256 of the generated structure -- the lever changed 158 of 1067 -- coordination.intact and elapsed_s, plus the sweep of record's sha256; and where each arm came from) |
| `v0418_sel_ab_off_bucket_report_honest.json.gz` | 53271 | `a03cac84d1a7666a` | tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py   (per arm: the honest bucket and the neutral ruler's verdict for every one of the 1146; ab_off is the v0.4.18 sweep's own rows. With v0.4.18-census/v0418_census_attribution_table.tsv.gz these re-derive the headline: tools/v0418/freeze_stage_selection.py --verify <this dir>) |
| `v0418_sel_ab_off_g_verdict.jsonl.gz` | 41364 | `dc569bec091d15ef` | tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py   (per arm: the honest bucket and the neutral ruler's verdict for every one of the 1146; ab_off is the v0.4.18 sweep's own rows. With v0.4.18-census/v0418_census_attribution_table.tsv.gz these re-derive the headline: tools/v0418/freeze_stage_selection.py --verify <this dir>) |
| `v0418_sel_ab_retarget_bucket_report_honest.json.gz` | 52587 | `b71075be4a0121f5` | tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py   (per arm: the honest bucket and the neutral ruler's verdict for every one of the 1146; ab_off is the v0.4.18 sweep's own rows. With v0.4.18-census/v0418_census_attribution_table.tsv.gz these re-derive the headline: tools/v0418/freeze_stage_selection.py --verify <this dir>) |
| `v0418_sel_ab_retarget_g_verdict.jsonl.gz` | 41808 | `297319d860b9053b` | tools/roundtrip_bucket_report.py --score honest ; tools/census/g_vs_input.py   (per arm: the honest bucket and the neutral ruler's verdict for every one of the 1146; ab_off is the v0.4.18 sweep's own rows. With v0.4.18-census/v0418_census_attribution_table.tsv.gz these re-derive the headline: tools/v0418/freeze_stage_selection.py --verify <this dir>) |
| `v0418_sel_ab_rows.tsv.gz` | 69428 | `e2114632634210fe` | tools/v0418/freeze_stage_selection.py   (MADE at staging time: per molecule and arm the sha256 of the generated structure -- the lever changed 158 of 1067 -- coordination.intact and elapsed_s, plus the sweep of record's sha256; and where each arm came from) |
| `v0418_sel_eta_ab_report_retarget.json` | 518 | `f1dfc630e635a35f` | tools/v0418/run_selection_ab.sh -> post_eta_ab.sh retarget -> eta_ab_report.py --on-arm retarget --base 4305,3883   (THE HARNESS A/B of OIN_ETA_RETARGET over all 1146 eta-bound molecules, commit 2f0b68d5: self-consistent +49/-4, VERIFIED +43/-5 -> 87.00% / 78.42% projected; excluding rows on the 300 s budget boundary +44/-4 and +40/-5. THE OFF ARM WAS BUILT from the v0.4.18 sweep of record, not run, so the report's noise-floor line is circular here. NO DEFAULT CHANGED: the lever ships OFF) |
| `v0418_sel_eta_ab_report_retarget.txt` | 3886 | `625512fa258c9e1e` | tools/v0418/run_selection_ab.sh -> post_eta_ab.sh retarget -> eta_ab_report.py --on-arm retarget --base 4305,3883   (THE HARNESS A/B of OIN_ETA_RETARGET over all 1146 eta-bound molecules, commit 2f0b68d5: self-consistent +49/-4, VERIFIED +43/-5 -> 87.00% / 78.42% projected; excluding rows on the 300 s budget boundary +44/-4 and +40/-5. THE OFF ARM WAS BUILT from the v0.4.18 sweep of record, not run, so the report's noise-floor line is circular here. NO DEFAULT CHANGED: the lever ships OFF) |
| `v0418_sel_probe5_classes.json` | 1790 | `349f7a82a354680c` | tools/v0418/eta_path_probe.py --arms dg,retarget   (fresh-process SAMPLE: 63 oracle-headroom molecules + 40 verified controls. Honest round trip 22 -> 52 (+31 / -1), controls 40/40; lever OFF reproduces the v0.4.18 sweep on 103/103 structures. A SAMPLE, NOT A VERDICT) |
| `v0418_sel_retarget_probe.txt` | 10307 | `8f83c5093469e7ab` | tools/v0418/eta_path_probe.py --arms dg,retarget   (fresh-process SAMPLE: 63 oracle-headroom molecules + 40 verified controls. Honest round trip 22 -> 52 (+31 / -1), controls 40/40; lever OFF reproduces the v0.4.18 sweep on 103/103 structures. A SAMPLE, NOT A VERDICT) |
| `v0418_sel_selection_sim.json` | 597 | `dbd527e750b01a9c` | tools/v0418/selection_sim.py   (OFFLINE BOUND, and it CAN express a loss: the L2 A/B left three structures per eta molecule with KNOWN verdicts. An input-free 'every binding atom the string DECLARES is inside the contact cutoff' choice between them = 668 -> 723 verified (+57 / -2), oracle 740. Adding SURPLUS contacts as a criterion costs 20 losses) |
| `v0418_sel_selection_sim.txt` | 1413 | `de2fb3d1539da24a` | tools/v0418/selection_sim.py   (OFFLINE BOUND, and it CAN express a loss: the L2 A/B left three structures per eta molecule with KNOWN verdicts. An input-free 'every binding atom the string DECLARES is inside the contact cutoff' choice between them = 668 -> 723 verified (+57 / -2), oracle 740. Adding SURPLUS contacts as a criterion costs 20 losses) |

Source paths at harvest time:

```
v0418_sel_ab_commits.tsv  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_ab_commits.tsv
v0418_sel_ab_off_bucket_report_honest.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_ab_off_bucket_report_honest.json.gz
v0418_sel_ab_off_g_verdict.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_ab_off_g_verdict.jsonl.gz
v0418_sel_ab_retarget_bucket_report_honest.json.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_ab_retarget_bucket_report_honest.json.gz
v0418_sel_ab_retarget_g_verdict.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_ab_retarget_g_verdict.jsonl.gz
v0418_sel_ab_rows.tsv.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_ab_rows.tsv.gz
v0418_sel_eta_ab_report_retarget.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_eta_ab_report_retarget.json
v0418_sel_eta_ab_report_retarget.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_eta_ab_report_retarget.txt
v0418_sel_probe5_classes.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_probe5_classes.json
v0418_sel_retarget_probe.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_retarget_probe.txt
v0418_sel_selection_sim.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_selection_sim.json
v0418_sel_selection_sim.txt  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.18-eta-selection/freeze/v0418_sel_selection_sim.txt
```
