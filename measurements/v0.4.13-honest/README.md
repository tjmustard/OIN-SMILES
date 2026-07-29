# `measurements/v0.4.13-honest` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `mirror_armA_promoted.json` | 20307 | `34169b792696f63c` | tools/mirror_audit_donor_fold.py --dataset <cohort-v0.4.5-5k> --n 250 --seed 7   (arm A: promoted; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- mixed cat+photo draw, reads 33 collapses -> 0 |
| `mirror_armB_noveto.json` | 20538 | `6cd6ca2ce98df783` | tools/mirror_audit_donor_fold.py --dataset <cohort-v0.4.5-5k> --n 250 --seed 7   (arm B: noveto; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- mixed cat+photo draw, reads 33 collapses -> 0 |
| `mirror_cat_noveto.json` | 20678 | `2b3284831ad3d262` | tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 7   (noveto; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- CAT-ONLY draw, reproduces v0.4.12's published 19 -> 0 with achiral unmoved at 157 |
| `mirror_cat_promoted.json` | 20530 | `83ff2c17de5d0aee` | tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 7   (promoted; the noveto arm sets OIN_FOLD_PARITY_VETO=0) -- CAT-ONLY draw, reproduces v0.4.12's published 19 -> 0 with achiral unmoved at 157 |
| `prefilter_AROHIA_comp_0.json` | 1219 | `e0d6b30d0a9091b8` | tools/prefilter_prevalence.py --xyz <input> | --cohort <dir>   (OIN_PREFILTER_ADVISORY two-arm; needs a QUIET machine — it reports a latency cost) |

Source paths at harvest time:

```
mirror_armA_promoted.json  <-  <SCRATCH>
mirror_armB_noveto.json  <-  <SCRATCH>
mirror_cat_noveto.json  <-  <SCRATCH>
mirror_cat_promoted.json  <-  <SCRATCH>
prefilter_AROHIA_comp_0.json  <-  <SCRATCH>
```
