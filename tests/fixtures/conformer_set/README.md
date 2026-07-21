# Conformer-invariance test set

Size-stratified sample of transition-metal complexes that round-trip cleanly
(XYZ -> OIN -> XYZ) in the sweep, used to verify that multiple conformers of a
structure collapse to the same canonical OIN-SMILES. The two Pt fixtures
(CisPlatin/TransPlatin) are always included as small, fast anchors and live in
`tests/fixtures/` (referenced in place, not copied here).

Regenerate deterministically with (point `--dataset-dir` at the gitignored
dataset in the main checkout):

```
python tools/select_conformer_test_set.py --n 30 --seed 42 --pool both --max-heavy 90 --max-per-metal 3
```

- pool: `both` (v0.4.0 quick success ∩ capstone-v042 accuracy-clean)
- seed: `42`  ·  max heavy atoms: `90`  ·  max per metal: `3`
- selected: **30** structures

## Heavy-atom strata (pool share vs. selected)

| bin | pool share | selected |
|---|---|---|
| <=20 | 7.8% | 4 |
| 21-30 | 30.6% | 8 |
| 31-40 | 29.4% | 8 |
| 41-50 | 19.1% | 5 |
| 51-75 | 12.5% | 4 |
| 76+ | 0.6% | 1 |

## Metal coverage (28 distinct)

Pt × 2, Ni × 2, Cu × 1, Au × 1, Pd × 1, Ti × 1, Hg × 1, Zn × 1, Rh × 1, Ag × 1, Fe × 1, Os × 1, Ru × 1, V × 1, Mo × 1, Cd × 1, Cr × 1, Co × 1, Ir × 1, W × 1, Mn × 1, Hf × 1, Zr × 1, Sc × 1, Nb × 1, Re × 1, Y × 1, Ta × 1

## Structures

| molecule | metal | heavy | total | mass | bin | charge | mult | source |
|---|---|---|---|---|---|---|---|---|
| ALIWOJ_comp_1 | Cu | 3 | 3 | 134.5 | <=20 | -1 | 1 | both |
| XIZFAQ_comp_0 | Au | 3 | 3 | 267.9 | <=20 | -1 | 1 | both |
| CisPlatin | Pt | 5 | 11 | 300.1 | <=20 | 0 | 1 | fixture |
| TransPlatin | Pt | 5 | 11 | 300.1 | <=20 | 0 | 1 | fixture |
| QEXRAP_comp_0 | Pd | 21 | 29 | 412.6 | 21-30 | 0 | 1 | both |
| MAZJED_comp_0 | Ti | 23 | 55 | 400.4 | 21-30 | 0 | 1 | both |
| OPAGES_comp_0 | Hg | 23 | 28 | 512.8 | 21-30 | 0 | 1 | both |
| USIVEZ_comp_0 | Zn | 25 | 37 | 373.7 | 21-30 | 0 | 1 | both |
| EGICON_comp_0 | Rh | 26 | 51 | 440.3 | 21-30 | 0 | 1 | both |
| WAGFOC_comp_0 | Ag | 28 | 36 | 568.0 | 21-30 | 0 | 1 | both |
| BEPCAC_comp_0 | Ni | 29 | 57 | 443.2 | 21-30 | 0 | 1 | both |
| QUSJAQ_comp_0 | Fe | 29 | 53 | 468.4 | 21-30 | 0 | 1 | both |
| XIXPEB_comp_0 | Os | 31 | 38 | 730.3 | 31-40 | -1 | 1 | both |
| RUKMEP_comp_0 | Ru | 32 | 52 | 523.6 | 31-40 | 1 | 1 | both |
| KARQUQ_comp_0 | V | 34 | 71 | 504.5 | 31-40 | 0 | 1 | both |
| FUXMUH_comp_0 | Mo | 35 | 47 | 720.6 | 31-40 | 0 | 1 | both |
| XIRLOA_comp_0 | Cd | 35 | 53 | 613.7 | 31-40 | 0 | 1 | both |
| NODNOL_comp_0 | Cr | 37 | 67 | 546.6 | 31-40 | 0 | 1 | both |
| PAGLUH_comp_0 | Co | 37 | 52 | 594.7 | 31-40 | 0 | 1 | both |
| COKGAN_comp_0 | Ir | 38 | 57 | 690.7 | 31-40 | 0 | 1 | both |
| AFIZEV_comp_0 | W | 41 | 60 | 831.1 | 41-50 | 0 | 1 | both |
| ZIGCIE_comp_0 | Mn | 43 | 72 | 642.5 | 41-50 | 0 | 1 | both |
| CETDAI_comp_0 | Hf | 47 | 112 | 814.5 | 41-50 | 0 | 1 | both |
| SINCIE_comp_0 | Zr | 49 | 81 | 708.0 | 41-50 | 0 | 1 | both |
| UKAMOJ_comp_0 | Sc | 49 | 73 | 675.6 | 41-50 | 0 | 1 | both |
| XEXMEV_comp_0 | Nb | 51 | 98 | 800.5 | 51-75 | 0 | 1 | both |
| NOVVUR_comp_0 | Re | 55 | 63 | 964.4 | 51-75 | 0 | 1 | both |
| YEYBEK_comp_0 | Y | 58 | 121 | 849.0 | 51-75 | 0 | 1 | both |
| CETPOI_comp_0 | Ta | 59 | 116 | 994.7 | 51-75 | 0 | 1 | both |
| YIJJEH_comp_0 | Ni | 77 | 147 | 1058.0 | 76+ | 0 | 1 | both |
