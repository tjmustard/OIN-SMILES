# `measurements/v0.4.17-census` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `attribution_summary.json` | 4450 | `77d0eb273709ff2a` | tools/census/attribution_table.py --gz   (census C4: the join -- one fault per molecule, first rule that fires, passes included; bucket column reproduces bucket_report_honest; UNATTRIBUTED printed as a number) |
| `attribution_table.tsv.gz` | 96342 | `9c3c3466da45c6b7` | tools/census/attribution_table.py --gz   (census C4: the join -- one fault per molecule, first rule that fires, passes included; bucket column reproduces bucket_report_honest; UNATTRIBUTED printed as a number) |
| `census_e_selfconsistency.jsonl.gz` | 268969 | `7818bba80902bd36` | tools/census/e_selfconsistency.py / tools/census/veto_probe.py   (census C2: eleven encodes per input -- identity, rewrite, rotation, 3 renumberings, 3 x 0.02 A noise, mirror; the veto probe re-encodes with OIN_FOLD_PARITY_VETO and the fold off to name the mechanism behind achiral molecules with two strings) |
| `census_g_verdict.jsonl.gz` | 185038 | `7755d4c4aac87a17` | tools/census/g_vs_input.py [--control mirror] [--crosstab]   (census C1: the neutral ruler -- element-labelled distance graph + donor-direction spheres + signed volumes, NO oinsmiles import -- judges every generated structure against its input; the mirror control gives every input's ruler chirality) |
| `census_g_verdict_control-mirror.jsonl.gz` | 171626 | `6fb5e96856e040ed` | tools/census/g_vs_input.py [--control mirror] [--crosstab]   (census C1: the neutral ruler -- element-labelled distance graph + donor-direction spheres + signed volumes, NO oinsmiles import -- judges every generated structure against its input; the mirror control gives every input's ruler chirality) |
| `census_parseback.jsonl.gz` | 156990 | `74b1f55c8b937f99` | tools/census/string_sufficiency.py {parseback [--side gen] | collide | pflags --charge-probe | summary}   (census C3: smiles_1 read back to a graph with NO 3D vs the ruler's input graph; natural-twin collision scan; perception flags with the stated charge honoured through a patched entry point) |
| `census_parseback_gen.jsonl.gz` | 165697 | `45e06ac5e9bca85e` | tools/census/string_sufficiency.py {parseback [--side gen] | collide | pflags --charge-probe | summary}   (census C3: smiles_1 read back to a graph with NO 3D vs the ruler's input graph; natural-twin collision scan; perception flags with the stated charge honoured through a patched entry point) |
| `census_pflags.jsonl.gz` | 68206 | `06ee726080a683f0` | tools/census/string_sufficiency.py {parseback [--side gen] | collide | pflags --charge-probe | summary}   (census C3: smiles_1 read back to a graph with NO 3D vs the ruler's input graph; natural-twin collision scan; perception flags with the stated charge honoured through a patched entry point) |
| `census_veto_probe.jsonl.gz` | 63221 | `923ca4484af3b41b` | tools/census/e_selfconsistency.py / tools/census/veto_probe.py   (census C2: eleven encodes per input -- identity, rewrite, rotation, 3 renumberings, 3 x 0.02 A noise, mirror; the veto probe re-encodes with OIN_FOLD_PARITY_VETO and the fold off to name the mechanism behind achiral molecules with two strings) |
| `census_veto_probe_chiral.jsonl.gz` | 105177 | `4339434793c09e78` | tools/census/e_selfconsistency.py / tools/census/veto_probe.py   (census C2: eleven encodes per input -- identity, rewrite, rotation, 3 renumberings, 3 x 0.02 A noise, mirror; the veto probe re-encodes with OIN_FOLD_PARITY_VETO and the fold off to name the mechanism behind achiral molecules with two strings) |
| `census_veto_probe_renumber.jsonl.gz` | 70725 | `e3cb3b2957ec3526` | tools/census/e_selfconsistency.py / tools/census/veto_probe.py   (census C2: eleven encodes per input -- identity, rewrite, rotation, 3 renumberings, 3 x 0.02 A noise, mirror; the veto probe re-encodes with OIN_FOLD_PARITY_VETO and the fold off to name the mechanism behind achiral molecules with two strings) |
| `collisions.json` | 41528 | `24d0373f5a9b98e3` | tools/census/string_sufficiency.py {parseback [--side gen] | collide | pflags --charge-probe | summary}   (census C3: smiles_1 read back to a graph with NO 3D vs the ruler's input graph; natural-twin collision scan; perception flags with the stated charge honoured through a patched entry point) |
| `g_crosstab.json` | 1953 | `43cb437327b138b6` | tools/census/g_vs_input.py [--control mirror] [--crosstab]   (census C1: the neutral ruler -- element-labelled distance graph + donor-direction spheres + signed volumes, NO oinsmiles import -- judges every generated structure against its input; the mirror control gives every input's ruler chirality) |
| `mirror_probe.json` | 35208 | `b95b5b6fde0d079c` | tools/census/mirror_probe.py   (census C1: does E(mirror x) differ from E(x) on a sample of ruler-achiral vs ruler-chiral inputs -- the estimate C2 then measured over the whole cohort) |
| `string_sufficiency_summary.json` | 8006 | `567166993cc27a17` | tools/census/string_sufficiency.py {parseback [--side gen] | collide | pflags --charge-probe | summary}   (census C3: smiles_1 read back to a graph with NO 3D vs the ruler's input graph; natural-twin collision scan; perception flags with the stated charge honoured through a patched entry point) |

Source paths at harvest time:

```
attribution_summary.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/attribution_summary.json
attribution_table.tsv.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/attribution_table.tsv.gz
census_e_selfconsistency.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_e_selfconsistency.jsonl.gz
census_g_verdict.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_g_verdict.jsonl.gz
census_g_verdict_control-mirror.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_g_verdict_control-mirror.jsonl.gz
census_parseback.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_parseback.jsonl.gz
census_parseback_gen.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_parseback_gen.jsonl.gz
census_pflags.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_pflags.jsonl.gz
census_veto_probe.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_veto_probe.jsonl.gz
census_veto_probe_chiral.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_veto_probe_chiral.jsonl.gz
census_veto_probe_renumber.jsonl.gz  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/census_veto_probe_renumber.jsonl.gz
collisions.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/collisions.json
g_crosstab.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/g_crosstab.json
mirror_probe.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/mirror_probe.json
string_sufficiency_summary.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-census/string_sufficiency_summary.json
```
