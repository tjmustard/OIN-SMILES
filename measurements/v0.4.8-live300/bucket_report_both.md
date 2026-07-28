# Round-trip bucket report -- SCORED vs HONEST, same molecules

Generated 2026-07-27T00:55:15 by tools/roundtrip_bucket_report.py --score both

Both columns classify the **same** reports with the **same** key and the **same** `status` gate.
The single variable is which round-trip string the verdict reads:

| column | string | what it is |
|---|---|---|
| scored | `smiles_2` | `get_oin_string(gen_result.mol, coords)` -- the generator's own bond graph. Asserts bonds the coordinates do not support; drops stereo they do. |
| honest | `smiles_2_indep` | a full `XYZToSMILES().convert()` of the generated XYZ -- bonds *and* stereo re-derived from coordinates alone. |

## Buckets

| bucket | scored | % | honest | % | delta |
|---|---:|---:|---:|---:|---:|
| byte_exact | 246 | 82.00% | 213 | 71.00% | -33 |
| key_equal | 28 | 9.33% | 33 | 11.00% | +5 |
| facmer_divergent | 0 | 0.00% | 1 | 0.33% | +1 |
| structural | 1 | 0.33% | 27 | 9.00% | +26 |
| hard_fail | 25 | 8.33% | 26 | 8.67% | +1 |
| encode_fail | 0 | 0.00% | 0 | 0.00% | +0 |
| **total** | **300** | | **300** | | |

## Transition matrix -- where every molecule went

| scored bucket | honest bucket | n |
|---|---|---:|
| byte_exact | byte_exact | 212 |
| hard_fail | hard_fail | 25 |
| key_equal | key_equal | 24 |
| byte_exact | structural | 23  **moved** |
| byte_exact | key_equal | 9  **moved** |
| key_equal | structural | 3  **moved** |
| byte_exact | facmer_divergent | 1  **moved** |
| byte_exact | hard_fail | 1  **moved** |
| key_equal | byte_exact | 1  **moved** |
| structural | structural | 1 |

## Molecules whose bucket moved (38)

- `BOCHIM_comp_0`  byte_exact -> structural
- `BODBOM_comp_0`  byte_exact -> facmer_divergent
- `CAZXIM_comp_0`  byte_exact -> structural
- `DILZOQ_comp_2`  byte_exact -> structural
- `DOKROM_comp_0`  byte_exact -> structural
- `EBUGEO_comp_0`  byte_exact -> key_equal
- `EWALET_comp_0`  byte_exact -> structural
- `FECSUC_comp_0`  byte_exact -> structural
- `GAMFOS_comp_0`  byte_exact -> structural
- `GANPIV_comp_0`  byte_exact -> key_equal
- `GIYBAU_comp_0`  byte_exact -> key_equal
- `HAVNUQ_comp_0`  byte_exact -> key_equal
- `HOBBUY_comp_0`  byte_exact -> hard_fail
- `ICUCAN_comp_0`  byte_exact -> structural
- `ICUYUC_comp_0`  byte_exact -> structural
- `JOWBOP_comp_0`  byte_exact -> structural
- `JUCCUH_comp_0`  byte_exact -> structural
- `KAFYUL_comp_0`  byte_exact -> key_equal
- `KIRVOY_comp_0`  byte_exact -> structural
- `KUJMUX_comp_0`  byte_exact -> structural
- `LABJEG_comp_0`  byte_exact -> structural
- `LASFER_comp_0`  byte_exact -> structural
- `LEYTIT_comp_0`  key_equal -> structural
- `MARXAH_comp_0`  key_equal -> byte_exact
- `MUPGEK_comp_0`  byte_exact -> key_equal
- `NOPXUM_comp_0`  byte_exact -> key_equal
- `OCIKAO_comp_0`  byte_exact -> structural
- `QOQMUF_comp_1`  byte_exact -> structural
- `QUHKIO_comp_0`  byte_exact -> structural
- `TIXWEC_comp_0`  byte_exact -> structural
- `UGINOO_comp_0`  key_equal -> structural
- `WIXLUL_comp_0`  byte_exact -> structural
- `XAXBAZ_comp_0`  byte_exact -> structural
- `XOQSON_comp_0`  byte_exact -> key_equal
- `XUZMIQ_comp_0`  byte_exact -> key_equal
- `YAVHIN_comp_1`  byte_exact -> structural
- `YIMVOG_comp_0`  key_equal -> structural
- `ZUFZUX_comp_0`  byte_exact -> structural

---

# v0.4.4 Round-Trip Bucket Report

Generated 2026-07-27T00:55:15 by tools/roundtrip_bucket_report.py
from 300 individual reports in `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.8-live300/individual_reports`,
classified with the v0.4.4 fac/mer-aware `oin.compare` key.

Round-trip string read from `smiles_2_indep` -- independent re-perception of the generated XYZ (**HONEST**).

## Buckets

| bucket | count | % |
|---|---:|---:|
| byte_exact | 213 | 71.00% |
| key_equal | 33 | 11.00% |
| facmer_divergent | 1 | 0.33% |
| structural | 27 | 9.00% |
| hard_fail | 26 | 8.67% |
| encode_fail | 0 | 0.00% |
| **total** | **300** | **100.00%** |

### key_equal sub-split (benign canonicalization reclaimed)

| subclass | count |
|---|---:|
| slot_renumber | 28 |
| rdkit_canonical | 5 |

## elapsed_s percentiles

| subset | n | p50 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|---:|
| overall | 300 | 7.5 | 123.6 | 300.2 | 300.5 | 448.8 |
| eta subset | 69 | 21.8 | 300.3 | 300.4 | 300.5 | 300.6 |

## fac/mer-divergent (newly-caught isomer errors -- SL2/SL3 target) (1)

- `BODBOM_comp_0`

## hard-fail (SL4 worklist) (26)

- `ATAGUZ_comp_2`
- `CAKBEW_comp_0`
- `DAYCUE_comp_0`
- `DOFCAE_comp_0`
- `EBUFOX_comp_0`
- `ECIGAZ_comp_0`
- `ESARUM_comp_0`
- `GADKUT_comp_0`
- `GURKUA_comp_0`
- `GUXMAP_comp_0`
- `HOBBUY_comp_0`
- `JINMID_comp_0`
- `LIMVEL_comp_0`
- `LIRFOJ_comp_0`
- `MAZKII_comp_0`
- `NAGQOE_comp_0`
- `NEFNER_comp_0`
- `NIHHAP_comp_0`
- `PDTPOR_comp_0`
- `RAJNIC_comp_0`
- `REKFAQ_comp_0`
- `RULMOA_comp_2`
- `SOQKEQ_comp_0`
- `TOGKOR_comp_0`
- `USAPOU_comp_0`
- `ZEKKIL_comp_0`

## encode-fail (SL5 worklist) (0)


## structural mismatch (27)

- `BOCHIM_comp_0`
- `BUWSUI_comp_0`
- `CAZXIM_comp_0`
- `DILZOQ_comp_2`
- `DOKROM_comp_0`
- `EWALET_comp_0`
- `FECSUC_comp_0`
- `GAMFOS_comp_0`
- `ICUCAN_comp_0`
- `ICUYUC_comp_0`
- `JOWBOP_comp_0`
- `JUCCUH_comp_0`
- `KIRVOY_comp_0`
- `KUJMUX_comp_0`
- `LABJEG_comp_0`
- `LASFER_comp_0`
- `LEYTIT_comp_0`
- `OCIKAO_comp_0`
- `QOQMUF_comp_1`
- `QUHKIO_comp_0`
- `TIXWEC_comp_0`
- `UGINOO_comp_0`
- `WIXLUL_comp_0`
- `XAXBAZ_comp_0`
- `YAVHIN_comp_1`
- `YIMVOG_comp_0`
- `ZUFZUX_comp_0`
