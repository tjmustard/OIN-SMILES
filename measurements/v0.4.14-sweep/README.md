# `measurements/v0.4.14-sweep` — frozen comparison artifacts

Written by `tools/harvest_measurements.py`. **Do not hand-edit** — rerun the tool.

| file | bytes | sha256 | produced by |
|---|---:|---|---|
| `RUN.md` | 2224 | `faf9bf13298df827` | hand-written provenance for the v0.4.14 baseline sweep: 6 shards 1-BASED, --mol-timeout 300, shipped lever defaults, BLAS threads capped to 1, and why systemd's OOMPolicy could not be used on a --scope |
| `attach_class_audit.json` | 126542 | `bd77288b962d4fda` | tools/attach_class_audit.py --results-dir <frozen sweep>   (MEDZUR/GAVSED split, with the byte_exact control arm) |
| `bucket_report.md` | 6536 | `2cc28b13f5bec6ec` | tools/roundtrip_bucket_report.py --results-dir <dir> |
| `bucket_report_PASS1_authoritative.md` | 15221 | `c6cda1d0f3ca5434` | tools/roundtrip_bucket_report.py --results-dir <sweep> --score honest, frozen   (THE v0.4.14 ABSOLUTE BASELINE: byte_exact 3858/5000 = 77.16%, gap 22.84) |
| `bucket_report_both.md` | 44728 | `ce862e886a9fda1e` | tools/roundtrip_bucket_report.py --results-dir <dir> |
| `bucket_report_honest.md` | 15221 | `c6cda1d0f3ca5434` | tools/roundtrip_bucket_report.py --results-dir <dir> |
| `per_molecule_extract.tsv` | 336931 | `71ce135907c2d3db` | derived from results-v0.4.14-sweep/individual_reports (N=5000): molecule, status, tier_passed, bucket under BOTH scored and honest, subclass, and the NESTED metrics.elapsed_s. Re-derives the authoritative table and the runtime percentiles without the 261 MB run |
| `veto_outcomes.json` | 53414 | `6081f1b84f7964b8` | tools/veto_outcome_audit.py --sweep <frozen sweep> --dataset <cat> --dataset <photo>   (which of fold_parity's FIVE outcomes each reverted molecule got: 222/222 vetoed_collapse, 0 no_evidence, over 393/393 movers) |
| `veto_residue_chirality.json` | 19776 | `9305c171a378bd4b` | tools/veto_residue_chirality.py --outcomes veto_outcomes.json --sweep <frozen sweep>   (183/222 MIRROR_MATCH: the round trip built the ENANTIOMER) |

Source paths at harvest time:

```
RUN.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/RUN.md
attach_class_audit.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/attach_class_audit.json
bucket_report.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/bucket_report.md
bucket_report_PASS1_authoritative.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/bucket_report_PASS1_authoritative.md
bucket_report_both.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/bucket_report_both.md
bucket_report_honest.md  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/bucket_report_honest.md
per_molecule_extract.tsv  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/per_molecule_extract.tsv
veto_outcomes.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/veto_outcomes.json
veto_residue_chirality.json  <-  tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/veto_residue_chirality.json
```
