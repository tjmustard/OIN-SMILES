# SL4 hard-fail triage — FF-only re-run of the BASELINE cohort

Baseline `hard_fail` cohort: **306**  
Re-run reports found: **306**  (not re-run: 0)  

- **Reclaimed (now key-exact success): 51  (16.7% of cohort)**
- Still failed: 255

## Reclaim reason

| reason | count | eta | non-eta |
|---|---:|---:|---:|
| rmsd_demote | 28 | 23 | 5 |
| passed_outright | 23 | 14 | 9 |

## Still-failed populations (disjoint)

| population | count | eta | non-eta |
|---|---:|---:|---:|
| atom_count | 113 | 24 | 89 |
| no_conformer | 109 | 35 | 74 |
| gen_exception | 25 | 6 | 19 |
| string_mismatch | 5 | 5 | 0 |
| timeout | 3 | 0 | 3 |

## Reclaimed structure quality (diagnostic)

Clash-free (clash_count == 0): **46/51** (90.2%); mean clash_count 0.37.

## Eta split (SL2 overlap)

- Reclaimed: 37 eta / 14 non-eta
- Still failed: 70 eta / 185 non-eta (the eta residual is SL2's target; re-measure after SL2 lands)
