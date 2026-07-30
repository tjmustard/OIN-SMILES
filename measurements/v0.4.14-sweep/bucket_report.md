# v0.4.4 Round-Trip Bucket Report

Generated 2026-07-29T18:43:43 by tools/roundtrip_bucket_report.py
from 5000 individual reports in `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.14-sweep/individual_reports`,
classified with the v0.4.4 fac/mer-aware `oin.compare` key.

Round-trip string read from `smiles_2` -- `get_oin_string(gen_result.mol, coords)`, the generator's own bond graph (**scored**, historical).

## Buckets

| bucket | count | % |
|---|---:|---:|
| byte_exact | 4344 | 86.88% |
| key_equal | 363 | 7.26% |
| facmer_divergent | 1 | 0.02% |
| structural | 18 | 0.36% |
| hard_fail | 262 | 5.24% |
| encode_fail | 12 | 0.24% |
| **total** | **5000** | **100.00%** |

### key_equal sub-split (benign canonicalization reclaimed)

| subclass | count |
|---|---:|
| slot_renumber | 301 |
| rdkit_canonical | 62 |

## elapsed_s percentiles

| subset | n | p50 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|---:|
| overall | 5000 | 4.0 | 52.1 | 193.4 | 354.5 | 728.8 |
| eta subset | 1146 | 11.0 | 139.6 | 300.2 | 411.7 | 648.8 |

## fac/mer-divergent (newly-caught isomer errors -- SL2/SL3 target) (1)

- `TULTAX_comp_0`

## hard-fail (SL4 worklist) (262)

- `AKUMAV_comp_0`
- `ALEMOT_comp_0`
- `ALITEU_comp_0`
- `AMUKEZ_comp_0`
- `ATAGUZ_comp_0`
- `ATAGUZ_comp_2`
- `BAQFIJ_comp_0`
- `BAXWUS_comp_0`
- `BEKLUA_comp_0`
- `BEWSUT_comp_0`
- `BOLDOW_comp_0`
- `BUKCUG_comp_0`
- `BUWHAD_comp_1`
- `CAJZOD_comp_0`
- `CAKBEW_comp_0`
- `CEGBAU_comp_0`
- `CEMTIZ_comp_0`
- `CEQWAY_comp_0`
- `CEQWUS_comp_0`
- `CETFAJ_comp_0`
- `CETYAC_comp_0`
- `COFMUH_comp_0`
- `COGKEQ_comp_0`
- `COJJAP_comp_0`
- `COJJOD_comp_0`
- `CUCCIO_comp_0`
- `DAGNAD_comp_0`
- `DAJREN_comp_0`
- `DAMVAQ_comp_0`
- `DAQDAB_comp_0`
- `DAQDIJ_comp_0`
- `DAYCUE_comp_0`
- `DEJHEF_comp_0`
- `DIFPAM_comp_0`
- `DIRFIW_comp_0`
- `DOCPAO_comp_0`
- `DOFCAE_comp_0`
- `DOFCAE_comp_2`
- `EBAFUL_comp_0`
- `EBAGAS_comp_0`
- `EBUBAH_comp_0`
- `ECIGAZ_comp_0`
- `EQEROI_comp_0`
- `EWAMAQ_comp_2`
- `FARFAI_comp_0`
- `FEKKAJ_comp_0`
- `FIQPIG_comp_0`
- `FIYBIB_comp_0`
- `FOJJUM_comp_0`
- `FOKDIV_comp_0`
- `FOSNEI_comp_0`
- `GANFEI_comp_0`
- `GEKNUG_comp_0`
- `GIFVUN_comp_0`
- `GOBYAA_comp_0`
- `GORNUY_comp_0`
- `GURJUA_comp_0`
- `GURKOU_comp_0`
- `GURKUA_comp_0`
- `GUXMAP_comp_0`
- `HACDOG_comp_0`
- `HAPCUY_comp_0`
- `HICLAG_comp_0`
- `HIMFAM_comp_0`
- `HIMQOL_comp_0`
- `HITHAW_comp_0`
- `HOHTEF_comp_0`
- `HOQFUP_comp_0`
- `HOZZUT_comp_0`
- `HUTCIK_comp_0`
- `HUTCOQ_comp_0`
- `HUTMEO_comp_0`
- `ICEZUO_comp_0`
- `ICOXIK_comp_0`
- `IDICEF_comp_0`
- `IFAPUD_comp_0`
- `IJIXEG_comp_0`
- `ILELOA_comp_0`
- `INALIU_comp_0`
- `IRINIH_comp_0`
- `ITOQIS_comp_0`
- `IXAHIZ_comp_0`
- `JEFHAE_comp_0`
- `JESFIX_comp_0`
- `JINMID_comp_0`
- `JIRFEW_comp_0`
- `JIVRIP_comp_0`
- `JIVVIT_comp_0`
- `JUXPID_comp_0`
- `KADYAS_comp_0`
- `KAXPAA_comp_0`
- `KECJUA_comp_0`
- `KEYVIU_comp_0`
- `KIKRUU_comp_0`
- `KIKSAB_comp_0`
- `KIKWOQ_comp_0`
- `KIVGIH_comp_0`
- `KODHIU_comp_0`
- `LAMTAX_comp_0`
- `LEKCIP_comp_0`
- `LETXOY_comp_0`
- `LEXDUP_comp_0`
- `LIMVEL_comp_0`
- `LISNUX_comp_0`
- `LOJGEW_comp_0`
- `LOLROW_comp_0`
- `MAFRUJ_comp_0`
- `MAHTOE_comp_0`
- `MAZKII_comp_0`
- `MIBFEL_comp_0`
- `MITZAS_comp_0`
- `MOCHUH_comp_0`
- `MOSLEL_comp_0`
- `MUKGUW_comp_0`
- `MUXKAT_comp_0`
- `NAHNOB_comp_0`
- `NARYOV_comp_0`
- `NASZOY_comp_0`
- `NAYJIG_comp_0`
- `NEFNER_comp_0`
- `NEGVOL_comp_0`
- `NEXTIT_comp_0`
- `NIXFAE_comp_0`
- `NIXFIM_comp_0`
- `NODKUP_comp_0`
- `NODLAW_comp_0`
- `NODLEA_comp_0`
- `NOJWAN_comp_0`
- `NOXREZ_comp_0`
- `NOYTUS_comp_0`
- `NUTHER_comp_0`
- `OBABAX_comp_0`
- `ODUBAS_comp_0`
- `ODUJEC_comp_0`
- `OFADAA_comp_0`
- `ONOXOG_comp_0`
- `OPUCUZ_comp_0`
- `OQAPED_comp_0`
- `OQIFEA_comp_0`
- `OTOLAL_comp_0`
- `OWEWAQ_comp_0`
- `OWODUA_comp_0`
- `OZADAW_comp_0`
- `PAFZUS_comp_0`
- `PAQCAM_comp_0`
- `PAWJED_comp_0`
- `PDTPOR_comp_0`
- `PEKRAB_comp_0`
- `PICVII_comp_0`
- `PICVUT_comp_0`
- `PICWAA_comp_0`
- `PIJPAB_comp_0`
- `PIJWEK_comp_0`
- `PODZEO_comp_0`
- `PUMKUD_comp_0`
- `PURROJ_comp_0`
- `QAMJAP_comp_0`
- `QEWJOS_comp_0`
- `QEXHEH_comp_0`
- `QIDHEQ_comp_0`
- `QIDKUL_comp_0`
- `QIFHIY_comp_0`
- `QIYYED_comp_0`
- `QIYYON_comp_0`
- `QOCHAT_comp_0`
- `QOFTUA_comp_0`
- `QOJKUW_comp_0`
- `QONNAK_comp_0`
- `RAJNIC_comp_0`
- `RATXOZ_comp_0`
- `RAXJEH_comp_0`
- `RAYCEC_comp_0`
- `REMVIQ_comp_0`
- `RIRYOJ_comp_0`
- `RIWKEO_comp_0`
- `RIXSAU_comp_1`
- `ROLYIB_comp_0`
- `ROMSER_comp_0`
- `RONNOA_comp_0`
- `RONQET_comp_0`
- `RONQOD_comp_0`
- `RUKWOL_comp_0`
- `RULMOA_comp_2`
- `SEMPEF_comp_0`
- `SEMPOP_comp_0`
- `SEMVUD_comp_0`
- `SEQVEP_comp_0`
- `SERCOI_comp_0`
- `SICNIC_comp_0`
- `SIGYAL_comp_0`
- `SOQKEQ_comp_0`
- `SUJDAC_comp_0`
- `SUMQAT_comp_0`
- `SUNROK_comp_0`
- `SUVREH_comp_0`
- `SUXJOM_comp_0`
- `SUXJUS_comp_0`
- `TANROQ_comp_0`
- `TEGVOR_comp_0`
- `TEGWOS_comp_0`
- `TEPGED_comp_0`
- `TEQHOM_comp_0`
- `TEYCOQ_comp_0`
- `TIKVOB_comp_1`
- `TIPDIG_comp_0`
- `TIZLEU_comp_0`
- `TOQXUU_comp_0`
- `TOXDIV_comp_0`
- `TUVCIX_comp_0`
- `ULOFEG_comp_0`
- `ULOQIX_comp_0`
- `ULOQOD_comp_0`
- `ULORAQ_comp_0`
- `UPABUK_comp_0`
- `UQUXAG_comp_0`
- `UROGAK_comp_0`
- `USAPOU_comp_0`
- `UTEMIR_comp_0`
- `UTOJEU_comp_0`
- `UZEVAY_comp_0`
- `VADCEI_comp_0`
- `VEJXOZ_comp_0`
- `VEQFUT_comp_0`
- `VIMTET_comp_0`
- `VIQQAN_comp_0`
- `VOPQOJ_comp_0`
- `VUXTOZ_comp_0`
- `WAHXOV_comp_1`
- `WAMWUE_comp_0`
- `WECYOS_comp_0`
- `WERLOV_comp_0`
- `WODCIC_comp_0`
- `WOGBOK_comp_0`
- `WOLRIA_comp_0`
- `WOTQAZ_comp_0`
- `WUMQAY_comp_0`
- `XAJBIW_comp_0`
- `XAKCAP_comp_0`
- `XENNIO_comp_0`
- `XENZAS_comp_0`
- `XOMKUG_comp_0`
- `XOSCIT_comp_0`
- `XUHXAA_comp_0`
- `YARYOI_comp_0`
- `YIDFEZ_comp_0`
- `YIQVOJ_comp_0`
- `YIVLAQ_comp_0`
- `YOLHOX_comp_0`
- `YOQLEY_comp_0`
- `YOQMAT_comp_0`
- `YOSYEM_comp_0`
- `YUVXOE_comp_0`
- `YUXZOI_comp_0`
- `ZAYNEU_comp_0`
- `ZEKKIL_comp_0`
- `ZENZAW_comp_1`
- `ZEPPAO_comp_0`
- `ZESCEI_comp_0`
- `ZIQLOD_comp_0`
- `ZISKUJ_comp_0`
- `ZOCYIC_comp_0`
- `ZOLXOP_comp_0`

## encode-fail (SL5 worklist) (12)

- `BENVOG_comp_0`
- `GABSOT_comp_0`
- `KEMSUS_comp_0`
- `KESWUB_comp_0`
- `LUJJUX_comp_0`
- `NOGWOX_comp_0`
- `OMEYUA_comp_0`
- `SUGQOA_comp_0`
- `UROFUD_comp_0`
- `UVIWIG_comp_0`
- `XAWSET_comp_0`
- `ZEHQOU_comp_0`

## structural mismatch (18)

- `AWORUZ_comp_0`
- `BUWSUI_comp_0`
- `DAQLOY_comp_0`
- `GUYQUO_comp_0`
- `KIHHUG_comp_0`
- `LAGJIP_comp_0`
- `NEBVAU_comp_0`
- `NOEPOR_comp_0`
- `PEDPEW_comp_0`
- `PEDPOG_comp_0`
- `PORKUF_comp_0`
- `PUVWEK_comp_0`
- `QIWNIU_comp_0`
- `SEBPUL_comp_1`
- `TOCZES_comp_0`
- `VUPQEE_comp_0`
- `ZIPFOV_comp_0`
- `ZOYYUJ_comp_0`
