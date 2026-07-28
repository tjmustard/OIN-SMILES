# `measurements/v0.4.13-honest` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `FROZEN.md` | 5245 | `9e5374f248ec361d` | hand-written provenance record for a frozen sweep |
| `SOURCE` | 92 | `26b5634f295c8cd3` | hand-written provenance record for a frozen sweep |
| `attach_class_audit.json` | 126547 | `80caa8e88c6e919e` | tools/attach_class_audit.py --results-dir <frozen sweep>   (MEDZUR/GAVSED split, with the byte_exact control arm) |
| `bucket_report_PASS1_authoritative.md` | 894 | `e9c24f05972df55c` | tools/roundtrip_bucket_report.py --results-dir <dir> |
| `fold_key_invariance.json` | 6771 | `10ab0e0b11f18f63` | tools/fold_key_invariance.py --sweep <frozen sweep>   (does the fold ever change the round-trip KEY? 0 => generator-neutral => an offline re-score is exact and no sweep is owed) |
| `fold_transition_veto.json` | 21149 | `e12ae96f74eed622` | tools/fold_transition_sim.py --sweep <frozen sweep> --arm veto   --dataset <cat> --dataset <photo>   (ABSOLUTE roots: the relative default silently excludes every mover from a worktree) |
| `mirror_armA_promoted.json` | 20307 | `34169b792696f63c` | tools/mirror_audit_donor_fold.py --dataset <cohort> --n 250 --seed 7   (arm A: promoted; arm B sets OIN_FOLD_PARITY_VETO=0) |
| `mirror_armB_noveto.json` | 20538 | `6cd6ca2ce98df783` | tools/mirror_audit_donor_fold.py --dataset <cohort> --n 250 --seed 7   (arm B: noveto; arm B sets OIN_FOLD_PARITY_VETO=0) |

Source paths at harvest time:

```
FROZEN.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/FROZEN.md
SOURCE  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/SOURCE
attach_class_audit.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/attach_class_audit.json
bucket_report_PASS1_authoritative.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/bucket_report_PASS1_authoritative.md
fold_key_invariance.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/fold_key_invariance.json
fold_transition_veto.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/fold_transition_veto.json
mirror_armA_promoted.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/mirror_armA_promoted.json
mirror_armB_noveto.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/mirror_armB_noveto.json
```
