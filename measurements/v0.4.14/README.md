# `measurements/v0.4.14` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `atrisk_generator_ab.json` | 6751 | `94b4eb9f1b849f31` | tools/generator_ab_honest.py --lever OIN_RESONANCE_DONOR_FOLD (sampling pass, superseded by reso_full_ab.json) |
| `generator_ab_honest.json` | 7167 | `8b49290034bf631b` | tools/generator_ab_honest.py --lever OIN_RESONANCE_DONOR_FOLD (sampling pass, superseded by reso_full_ab.json) |
| `hekfel_honest_ab.json` | 994 | `08fde66e79675261` | tools/generator_ab_honest.py --lever OIN_RESONANCE_DONOR_FOLD (sampling pass, superseded by reso_full_ab.json) |
| `mirror_cat_control_resoff.json` | 20530 | `83ff2c17de5d0aee` | tools/mirror_audit_donor_fold.py --dataset <cat cohort> --n cat   (OIN_RESONANCE_DONOR_FOLD=control_resoff; the `movers` draw is the MOVER-ENRICHED cohort at 179/179 coverage -- the `cat` draw holds only 1 mover in 250 = 0.4% and is a general-population control, NOT evidence about this lever) |
| `mirror_cat_resonance.json` | 20530 | `83ff2c17de5d0aee` | tools/mirror_audit_donor_fold.py --dataset <cat cohort> --n cat   (OIN_RESONANCE_DONOR_FOLD=resonance; the `movers` draw is the MOVER-ENRICHED cohort at 179/179 coverage -- the `cat` draw holds only 1 mover in 250 = 0.4% and is a general-population control, NOT evidence about this lever) |
| `mirror_movers_control_resoff.json` | 14343 | `0745893b87b3110c` | tools/mirror_audit_donor_fold.py --dataset <movers cohort> --n movers   (OIN_RESONANCE_DONOR_FOLD=control_resoff; the `movers` draw is the MOVER-ENRICHED cohort at 179/179 coverage -- the `cat` draw holds only 1 mover in 250 = 0.4% and is a general-population control, NOT evidence about this lever) |
| `mirror_movers_resonance.json` | 14343 | `0745893b87b3110c` | tools/mirror_audit_donor_fold.py --dataset <movers cohort> --n movers   (OIN_RESONANCE_DONOR_FOLD=resonance; the `movers` draw is the MOVER-ENRICHED cohort at 179/179 coverage -- the `cat` draw holds only 1 mover in 250 = 0.4% and is a general-population control, NOT evidence about this lever) |
| `prefilter_AROHIA_wiring_gate.json` | 1219 | `6aab6d67eef59f9b` | tools/prefilter_prevalence.py --xyz <input> | --cohort <dir>   (OIN_PREFILTER_ADVISORY two-arm; needs a QUIET machine — it reports a latency cost) |
| `prefilter_prevalence.json` | 31150 | `460677eae986737c` | tools/prefilter_prevalence.py --xyz <input> | --cohort <dir>   (OIN_PREFILTER_ADVISORY two-arm; needs a QUIET machine — it reports a latency cost) |
| `reso_full_ab.json` | 60304 | `1b4655756806fa02` | tools/generator_ab_honest.py --lever OIN_RESONANCE_DONOR_FOLD over all 182 affected   (THE v0.4.14 HEADLINE: 78 gains, 7 losses, net +71 = +1.42 pts, n=182 of 182) |
| `reso_movers_exact.json` | 2076 | `2b3c3c74752c9801` | tools/lever_string_movers.py --lever OIN_RESONANCE_DONOR_FOLD --holding OIN_CANONICAL_DONOR_FOLD --holding OIN_FOLD_PARITY_VETO   (93 of 5000 move encode(input) -- the coordinate-derived affected population) |
| `resonance_key_invariance.json` | 6902 | `3e3ba28460b73be4` | tools/fold_key_invariance.py --sweep <frozen sweep> --lever OIN_RESONANCE_DONOR_FOLD --holding OIN_CANONICAL_DONOR_FOLD   (9669 compared, 228 moved, 0 keys changed. ⚠ bounds ACCEPTANCE, not embedding) |
| `resonance_transition.json` | 24182 | `bb4455c7579fae32` | tools/resonance_transition_sim.py --sweep <frozen sweep> --baseline-byte-exact 3794   (OFFLINE re-score. Reported +78/0 losses; SUPERSEDED by reso_full_ab.json -- an offline re-score cannot express a loss) |
| `v0413_atrisk_ab.json` | 13195 | `67fce9fa08911851` | tools/generator_ab_honest.py --lever OIN_CANONICAL_DONOR_FOLD over 40 of v0.4.13's 197 at-risk molecules (seed 13)   (6 losses = 15% => ~30 over the population => v0.4.13's true net ~+2.82, NOT +3.42) |
| `veto_outcomes.json` | 86571 | `df9312f92d35d796` | tools/veto_outcome_audit.py --sweep <frozen sweep> --dataset <cat> --dataset <photo>   (which of fold_parity's FIVE outcomes each reverted molecule got: 222/222 vetoed_collapse, 0 no_evidence, over 393/393 movers) |
| `veto_residue_chirality.json` | 18176 | `f4a096265472eac1` | tools/veto_residue_chirality.py --outcomes veto_outcomes.json --sweep <frozen sweep>   (183/222 MIRROR_MATCH: the round trip built the ENANTIOMER) |

Source paths at harvest time:

```
atrisk_generator_ab.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/atrisk_generator_ab.json
generator_ab_honest.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/generator_ab_honest.json
hekfel_honest_ab.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/hekfel_honest_ab.json
mirror_cat_control_resoff.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/mirror_cat_control_resoff.json
mirror_cat_resonance.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/mirror_cat_resonance.json
mirror_movers_control_resoff.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/mirror_movers_control_resoff.json
mirror_movers_resonance.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/mirror_movers_resonance.json
prefilter_AROHIA_wiring_gate.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane2/prefilter_AROHIA_wiring_gate.json
prefilter_prevalence.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane2/prefilter_prevalence.json
reso_full_ab.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/reso_full_ab.json
reso_movers_exact.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/reso_movers_exact.json
resonance_key_invariance.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/resonance_key_invariance.json
resonance_transition.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/resonance_transition.json
v0413_atrisk_ab.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/v0413_atrisk_ab.json
veto_outcomes.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/veto_outcomes.json
veto_residue_chirality.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-lane1/veto_residue_chirality.json
```
