# v0.4.4 Round-Trip Bucket Report

Generated 2026-07-26T20:21:56 by tools/roundtrip_bucket_report.py
from 5000 individual reports in `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep/individual_reports`,
classified with the v0.4.4 fac/mer-aware `oin.compare` key.

## Buckets

| bucket | count | % |
|---|---:|---:|
| byte_exact | 4140 | 82.80% |
| key_equal | 520 | 10.40% |
| facmer_divergent | 5 | 0.10% |
| structural | 20 | 0.40% |
| hard_fail | 300 | 6.00% |
| encode_fail | 15 | 0.30% |
| **total** | **5000** | **100.00%** |

### key_equal sub-split (benign canonicalization reclaimed)

| subclass | count |
|---|---:|
| slot_renumber | 459 |
| rdkit_canonical | 61 |

## elapsed_s percentiles

| subset | n | p50 | p90 | p95 | p99 | max |
|---|---:|---:|---:|---:|---:|---:|
| overall | 5000 | 7.2 | 106.0 | 300.1 | 300.8 | 759.9 |
| eta subset | 1146 | 23.7 | 300.1 | 300.3 | 315.4 | 511.5 |

## fac/mer-divergent (newly-caught isomer errors -- SL2/SL3 target) (5)

- `QEWJOS_comp_0`
- `ROMSER_comp_0`
- `SEQVEP_comp_0`
- `TEQHOM_comp_0`
- `TULTAX_comp_0`

## hard-fail (SL4 worklist) (300)

- `ADEZOY_comp_0`
- `AKUMAV_comp_0`
- `ALEMOT_comp_0`
- `ALITEU_comp_0`
- `AMUKEZ_comp_0`
- `ATAGUZ_comp_0`
- `ATAGUZ_comp_2`
- `BAJBAS_comp_0`
- `BAQFIJ_comp_0`
- `BAXWUS_comp_0`
- `BEDLII_comp_0`
- `BEKLUA_comp_0`
- `BEWSUT_comp_0`
- `BOLDOW_comp_0`
- `BUKCUG_comp_0`
- `BUWHAD_comp_1`
- `CAHQOT_comp_0`
- `CAJZOD_comp_0`
- `CAKBEW_comp_0`
- `CEGBAU_comp_0`
- `CEMTIZ_comp_0`
- `CEQWAY_comp_0`
- `CEQWUS_comp_0`
- `CETFAJ_comp_0`
- `CETYAC_comp_0`
- `CINXAA_comp_0`
- `COFMUH_comp_0`
- `COGKEQ_comp_0`
- `COJJAP_comp_0`
- `COJJOD_comp_0`
- `CUCCIO_comp_0`
- `DAGNAD_comp_0`
- `DAGNIL_comp_0`
- `DAMVAQ_comp_0`
- `DAYCUE_comp_0`
- `DEJHEF_comp_0`
- `DIFPAM_comp_0`
- `DIRFIW_comp_0`
- `DOCPAO_comp_0`
- `DOFCAE_comp_0`
- `DOFCAE_comp_2`
- `DOTBIA_comp_0`
- `EBAFUL_comp_0`
- `EBAGAS_comp_0`
- `EBUBAH_comp_0`
- `ECIGAZ_comp_0`
- `EQEROI_comp_0`
- `ESARUM_comp_0`
- `ETILUP_comp_0`
- `EWAMAQ_comp_2`
- `FARFAI_comp_0`
- `FEKKAJ_comp_0`
- `FIQPIG_comp_0`
- `FIYBIB_comp_0`
- `FOQGEY_comp_0`
- `GADKUT_comp_0`
- `GANFEI_comp_0`
- `GEKNUG_comp_0`
- `GIFVUN_comp_0`
- `GOBYAA_comp_0`
- `GORNUY_comp_0`
- `GURJUA_comp_0`
- `GURKOU_comp_0`
- `GURKUA_comp_0`
- `GURWAT_comp_0`
- `GUXMAP_comp_0`
- `HACDOG_comp_0`
- `HAPCUY_comp_0`
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
- `IFADUQ_comp_0`
- `IFAPUD_comp_0`
- `IGITOI_comp_0`
- `IJIXEG_comp_0`
- `ILELOA_comp_0`
- `INALIU_comp_0`
- `IRINIH_comp_0`
- `IRUPOC_comp_0`
- `ITOQIS_comp_0`
- `IXAHIZ_comp_0`
- `JAYLIE_comp_0`
- `JEFHAE_comp_0`
- `JEHNES_comp_0`
- `JESFIX_comp_0`
- `JINMID_comp_0`
- `JIRFEW_comp_0`
- `JIVRIP_comp_0`
- `JIVVIT_comp_0`
- `KADYAS_comp_0`
- `KAQDEL_comp_0`
- `KAQFEN_comp_0`
- `KAXPAA_comp_0`
- `KECJUA_comp_0`
- `KEYVIU_comp_0`
- `KIKRUU_comp_0`
- `KIKSAB_comp_0`
- `KIKWOQ_comp_0`
- `KIVGIH_comp_0`
- `KODHIU_comp_0`
- `LACNOU_comp_0`
- `LAGJIP_comp_0`
- `LAHVOG_comp_0`
- `LAMTAX_comp_0`
- `LARQED_comp_0`
- `LEKCIP_comp_0`
- `LETXOY_comp_0`
- `LIMVEL_comp_0`
- `LIRFOJ_comp_0`
- `LISNUX_comp_0`
- `LOJGEW_comp_0`
- `LUMDON_comp_0`
- `LUSBEI_comp_0`
- `LUYYOT_comp_0`
- `LUZDAL_comp_0`
- `MAFRUJ_comp_0`
- `MAHTOE_comp_0`
- `MAZKII_comp_0`
- `MECXEY_comp_0`
- `MIBFEL_comp_0`
- `MIFJIV_comp_1`
- `MIHLIY_comp_0`
- `MITZAS_comp_0`
- `MOCHUH_comp_0`
- `MOPVEU_comp_0`
- `MUKGUW_comp_0`
- `MUXKAT_comp_0`
- `NAGQOE_comp_0`
- `NAHNOB_comp_0`
- `NASZOY_comp_0`
- `NAYJIG_comp_0`
- `NEFNER_comp_0`
- `NEGVOL_comp_0`
- `NEVCOI_comp_0`
- `NEVCUO_comp_0`
- `NEXTIT_comp_0`
- `NIHHAP_comp_0`
- `NIXFAE_comp_0`
- `NIXFIM_comp_0`
- `NIYJUA_comp_0`
- `NODKUP_comp_0`
- `NODLAW_comp_0`
- `NODLEA_comp_0`
- `NOEPOR_comp_0`
- `NOJWAN_comp_0`
- `NOYREA_comp_0`
- `NOYTUS_comp_0`
- `NURTEA_comp_0`
- `NUTHER_comp_0`
- `OBABAX_comp_0`
- `ODUBAS_comp_0`
- `ODUJEC_comp_0`
- `OPUCUZ_comp_0`
- `OQAPED_comp_0`
- `OTOLAL_comp_0`
- `OWEWAQ_comp_0`
- `OWODUA_comp_0`
- `OZADAW_comp_0`
- `OZUZEQ_comp_0`
- `PACMEL_comp_0`
- `PAFZUS_comp_0`
- `PAQCAM_comp_0`
- `PAXJUT_comp_0`
- `PDTPOR_comp_0`
- `PEDPEW_comp_0`
- `PEDPOG_comp_0`
- `PEKRAB_comp_0`
- `PICVII_comp_0`
- `PICVUT_comp_0`
- `PICWAA_comp_0`
- `PIJPAB_comp_0`
- `PIJWEK_comp_0`
- `PIKLOM_comp_0`
- `PODZEO_comp_0`
- `POZSIG_comp_0`
- `PUGVOC_comp_0`
- `PUMKUD_comp_0`
- `PUMLEP_comp_0`
- `PUQRAW_comp_0`
- `PURROJ_comp_0`
- `PUVWEK_comp_0`
- `QAMJAP_comp_0`
- `QEBKUG_comp_0`
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
- `RAYCEC_comp_0`
- `REMVIQ_comp_0`
- `REVMEL_comp_0`
- `RIRYOJ_comp_0`
- `RIVNUG_comp_0`
- `RIWKEO_comp_0`
- `RIXSAU_comp_1`
- `RONNOA_comp_0`
- `RONQET_comp_0`
- `RONQOD_comp_0`
- `RUKVAU_comp_0`
- `RUKWOL_comp_0`
- `RULMOA_comp_2`
- `SEMPEF_comp_0`
- `SEMPOP_comp_0`
- `SEMVUD_comp_0`
- `SERCOI_comp_0`
- `SICNIC_comp_0`
- `SIGYAL_comp_0`
- `SOQHEL_comp_0`
- `SOQKEQ_comp_0`
- `SUJDAC_comp_0`
- `SUMQAT_comp_0`
- `SUNROK_comp_0`
- `SUSGOC_comp_0`
- `SUVREH_comp_0`
- `SUXJOM_comp_0`
- `SUXJUS_comp_0`
- `TANROQ_comp_0`
- `TEPGED_comp_0`
- `TEYCOQ_comp_0`
- `TIKVOB_comp_1`
- `TIPDIG_comp_0`
- `TIZLEU_comp_0`
- `TOGKOR_comp_0`
- `TOQXUU_comp_0`
- `TOWWUX_comp_0`
- `TOXDIV_comp_0`
- `TUVCIX_comp_0`
- `UDARIC_comp_0`
- `ULOFEG_comp_0`
- `ULOQIX_comp_0`
- `ULOQOD_comp_0`
- `UPABUK_comp_0`
- `UQUXAG_comp_0`
- `USAPOU_comp_0`
- `UTEMIR_comp_0`
- `UTOJEU_comp_0`
- `UZEVAY_comp_0`
- `VADCEI_comp_0`
- `VAJVOS_comp_0`
- `VEJXOZ_comp_0`
- `VEMQIP_comp_0`
- `VIMTET_comp_0`
- `VIQQAN_comp_0`
- `VOPQOJ_comp_0`
- `VUXTOZ_comp_0`
- `WAHXOV_comp_1`
- `WAMWUE_comp_0`
- `WAVDED_comp_0`
- `WAWMUF_comp_0`
- `WECYOS_comp_0`
- `WERLOV_comp_0`
- `WIBTAD_comp_0`
- `WODCIC_comp_0`
- `WOFGUT_comp_0`
- `WOLRIA_comp_0`
- `WUMQAY_comp_0`
- `XAJBIW_comp_0`
- `XAKCAP_comp_0`
- `XANKAA_comp_0`
- `XATXOH_comp_0`
- `XENNIO_comp_0`
- `XENZAS_comp_0`
- `XILCOM_comp_0`
- `XIVBIN_comp_0`
- `XOSCIT_comp_0`
- `XOYCOE_comp_0`
- `YARYOI_comp_0`
- `YIDFEZ_comp_0`
- `YIQVOJ_comp_0`
- `YIVLAQ_comp_0`
- `YOLHOX_comp_0`
- `YOQMAT_comp_0`
- `YOSYEM_comp_0`
- `YUVXOE_comp_0`
- `ZAYNEU_comp_0`
- `ZENZAW_comp_1`
- `ZESCEI_comp_0`
- `ZIQLOD_comp_0`
- `ZISKUJ_comp_0`
- `ZOCYIC_comp_0`
- `ZODBUP_comp_0`
- `ZOLXOP_comp_0`
- `ZUCKOZ_comp_1`
- `ZUGCAH_comp_0`
- `ZUYVEU_comp_0`

## encode-fail (SL5 worklist) (15)

- `BENVOG_comp_0`
- `GABSOT_comp_0`
- `HICLAG_comp_0`
- `KEMSUS_comp_0`
- `KESWUB_comp_0`
- `LEXDUP_comp_0`
- `LUJJUX_comp_0`
- `NOGWOX_comp_0`
- `NOXREZ_comp_0`
- `OMEYUA_comp_0`
- `SUGQOA_comp_0`
- `UROFUD_comp_0`
- `UVIWIG_comp_0`
- `XAWSET_comp_0`
- `ZEHQOU_comp_0`

## structural mismatch (20)

- `AWORUZ_comp_0`
- `BUWSUI_comp_0`
- `FOSNEI_comp_0`
- `GUYQUO_comp_0`
- `KIHHUG_comp_0`
- `LOLROW_comp_0`
- `MOSLEL_comp_0`
- `NEBVAU_comp_0`
- `OQIFEA_comp_0`
- `PAWJED_comp_0`
- `PORKUF_comp_0`
- `QIWNIU_comp_0`
- `TOCZES_comp_0`
- `ULORAQ_comp_0`
- `VUPQEE_comp_0`
- `XOMKUG_comp_0`
- `YUXZOI_comp_0`
- `ZEKKIL_comp_0`
- `ZIPFOV_comp_0`
- `ZOYYUJ_comp_0`
