# v0.4.4 Round-Trip Bucket Report

Generated 2026-07-25T15:12:37 by tools/roundtrip_bucket_report.py
from 13 individual reports in `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-reencode/lane3-13cases-base/individual_reports`,
classified with the v0.4.4 fac/mer-aware `oin.compare` key.

## Buckets

| bucket | count | % |
|---|---:|---:|
| byte_exact | 0 | 0.00% |
| key_equal | 9 | 69.23% |
| facmer_divergent | 1 | 7.69% |
| structural | 3 | 23.08% |
| hard_fail | 0 | 0.00% |
| encode_fail | 0 | 0.00% |
| **total** | **13** | **100.00%** |

### key_equal sub-split (benign canonicalization reclaimed)

| subclass | count |
|---|---:|
| winding_star_drift | 6 |
| rdkit_canonical | 3 |

## elapsed_s percentiles

| subset | n | p50 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|---:|
| overall | 13 | 67.0 | 315.9 | 338.6 | 349.0 | 351.6 |
| eta subset | 13 | 67.0 | 315.9 | 338.6 | 349.0 | 351.6 |

## fac/mer-divergent (newly-caught isomer errors -- SL2/SL3 target) (1)

- `FAHCIC_comp_0`

## hard-fail (SL4 worklist) (0)


## encode-fail (SL5 worklist) (0)


## structural mismatch (3)

- `ABETIK_comp_0`
- `IFICAD_comp_0`
- `OVUBEO_comp_0`
