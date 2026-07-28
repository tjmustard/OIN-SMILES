# `measurements/v0.4.8-honest` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `FROZEN.md` | 1165 | `a077e857d9ef05ab` | hand-written provenance record for a frozen sweep |
| `SOURCE` | 153 | `cb3524c3b73c2194` | hand-written provenance record for a frozen sweep |
| `bucket_report_PASS1_authoritative.md` | 45063 | `0a8b0204c6a18707` | tools/roundtrip_bucket_report.py --results-dir <dir> |
| `bucket_report_both.md` | 45063 | `0a8b0204c6a18707` | tools/roundtrip_bucket_report.py --results-dir <dir> |
| `bucket_report_honest.md` | 15043 | `bc39002104d688c6` | tools/roundtrip_bucket_report.py --results-dir <dir> |

Source paths at harvest time:

```
FROZEN.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/FROZEN.md
SOURCE  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/SOURCE
bucket_report_PASS1_authoritative.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/bucket_report_PASS1_authoritative.md
bucket_report_both.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/bucket_report_both.md
bucket_report_honest.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-honest/bucket_report_honest.md
```
