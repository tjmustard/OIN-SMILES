# `measurements/v0.4.12` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `ab_eta_accept_realpop.json` | 20312 | `606a0af20a630821` | tools/ab_accept_scored.py --cohort <cohort> --lever <lever> --timeout 150 --hard-cap 240 |
| `ab_eta_accept_stale_cohort.json` | 14530 | `9d8e084785fcf7de` | tools/ab_accept_scored.py --cohort <cohort> --lever <lever> --timeout 150 --hard-cap 240 |
| `cohort_eta_fail.json` | 1243 | `693249803c660dc0` | cohort manifest, built from the frozen sweep |
| `cohort_pilot_stale.json` | 841 | `481fec24734a4ab2` | cohort manifest, built from the frozen sweep |
| `cohort_slow100.json` | 11731 | `7032acabc859d398` | cohort manifest, built from the frozen sweep |
| `mirror_audit_seed07_veto_off.json` | 20678 | `2b3284831ad3d262` | tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 07   (OIN_FOLD_PARITY_VETO=off) |
| `mirror_audit_seed07_veto_on.json` | 20530 | `83ff2c17de5d0aee` | tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 07   (OIN_FOLD_PARITY_VETO=on) |
| `mirror_audit_seed11_veto_off.json` | 20517 | `7eb17d282d01bc54` | tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 11   (OIN_FOLD_PARITY_VETO=off) |
| `mirror_audit_seed11_veto_on.json` | 20370 | `9ee0ca619486bf92` | tools/mirror_audit_donor_fold.py --dataset <cat> --n 250 --seed 11   (OIN_FOLD_PARITY_VETO=on) |
| `transition_fold.json` | 47469 | `ee8b673500c82ed2` | tools/fold_transition_sim.py --sweep <frozen sweep> --arm fold |
| `transition_veto.json` | 21105 | `e12ae96f74eed622` | tools/fold_transition_sim.py --sweep <frozen sweep> --arm veto |

Source paths at harvest time:

```
ab_eta_accept_realpop.json  <-  <SCRATCH>
ab_eta_accept_stale_cohort.json  <-  <SCRATCH>
cohort_eta_fail.json  <-  <SCRATCH>
cohort_pilot_stale.json  <-  <SCRATCH>
cohort_slow100.json  <-  <SCRATCH>
mirror_audit_seed07_veto_off.json  <-  <SCRATCH>
mirror_audit_seed07_veto_on.json  <-  <SCRATCH>
mirror_audit_seed11_veto_off.json  <-  <SCRATCH>
mirror_audit_seed11_veto_on.json  <-  <SCRATCH>
transition_fold.json  <-  <SCRATCH>
transition_veto.json  <-  <SCRATCH>
```
