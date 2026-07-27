```
--- ARM A: default (independent confirm ON) ---
  [A-default] YIYGAP_comp_0        1.83s pass=True clash=0 worst=0.84
  [A-default] NOMMOU_comp_0       39.89s pass=True clash=0 worst=0.7579
  [A-default] MEDZUR_comp_0        7.25s pass=True clash=0 worst=0.7577
  [A-default] WIWRIE_comp_0       74.55s pass=True clash=0 worst=0.7539
  [A-default] AROHIA_comp_0       70.53s pass=False clash=0 worst=0.8152
  [A-default] YIZHIY_comp_0        4.54s pass=True clash=0 worst=0.7588
  [A-default] FEXYOZ_comp_0       12.14s pass=True clash=0 worst=0.7621
  [A-default] XIQKOY_comp_0       306.5s pass=False clash=- worst=-
  [A-default] ODEWID_comp_0        1.45s pass=True clash=0 worst=0.7502
  [A-default] ZITSIE_comp_0      134.15s pass=True clash=0 worst=0.7526
  [A-default] QESRUE_comp_0        3.59s pass=True clash=0 worst=0.8206
  [A-default] HEJXIF_comp_0      194.59s pass=True clash=0 worst=0.7632
  [A-default] RATPEK_comp_0       68.95s pass=True clash=0 worst=0.7599
  [A-default] KAQDOV_comp_0      242.51s pass=True clash=0 worst=0.7504
  [A-default] GAVSED_comp_0    KILLED at hard cap 330s
  [A-default] GAVSED_comp_0           ?s pass=False clash=- worst=-
  [A-default] DEJHEF_comp_0        0.23s pass=False clash=- worst=-
  [A-default] XUPTAF_comp_0         5.6s pass=True clash=0 worst=0.8072
  [A-default] DAKGON_comp_0       13.87s pass=True clash=0 worst=0.8164
  [A-default] QIDKUL_comp_0    KILLED at hard cap 330s
  [A-default] QIDKUL_comp_0           ?s pass=False clash=- worst=-
  [A-default] LIYXEY_comp_0        6.82s pass=True clash=0 worst=0.7513
  [A-default] POVPIA_comp_0       13.89s pass=True clash=16 worst=0.4344
  [A-default] YENDUS_comp_0    KILLED at hard cap 330s
  [A-default] YENDUS_comp_0    KILLED at hard cap 330s
  [A-default] YENDUS_comp_0           ?s pass=False clash=- worst=-

--- ARM B: OIN_ACCEPT_SCORED=1 ---
  [B-scored] YIYGAP_comp_0        4.05s pass=True clash=0 worst=0.84
  [B-scored] NOMMOU_comp_0        4.08s pass=True clash=0 worst=0.7601
  [B-scored] MEDZUR_comp_0        1.83s pass=True clash=0 worst=0.7705
  [B-scored] WIWRIE_comp_0        6.34s pass=True clash=0 worst=0.7559
  [B-scored] AROHIA_comp_0       77.35s pass=False clash=0 worst=0.8152
  [B-scored] YIZHIY_comp_0        5.88s pass=True clash=0 worst=0.7588
  [B-scored] FEXYOZ_comp_0        4.38s pass=True clash=0 worst=0.7659
  [B-scored] ZITSIE_comp_0        5.46s pass=True clash=0 worst=0.7532
  [B-scored] ODEWID_comp_0        2.21s pass=True clash=0 worst=0.7502
  [B-scored] HEJXIF_comp_0       13.01s pass=True clash=0 worst=0.7562
  [B-scored] QESRUE_comp_0        5.63s pass=True clash=0 worst=0.8206
  [B-scored] KAQDOV_comp_0        7.48s pass=True clash=0 worst=0.7629
  [B-scored] RATPEK_comp_0       11.78s pass=True clash=1 worst=0.7461
  [B-scored] GAVSED_comp_0        5.41s pass=True clash=0 worst=0.7611
  [B-scored] XIQKOY_comp_0    KILLED at hard cap 330s
  [B-scored] XIQKOY_comp_0           ?s pass=False clash=- worst=-
  [B-scored] DEJHEF_comp_0         0.2s pass=False clash=- worst=-
  [B-scored] XUPTAF_comp_0        3.95s pass=True clash=0 worst=0.8072
  [B-scored] DAKGON_comp_0        3.55s pass=True clash=1 worst=0.7283
  [B-scored] POVPIA_comp_0        9.67s pass=True clash=0 worst=0.75
  [B-scored] QIDKUL_comp_0      165.95s pass=True clash=0 worst=0.7517
  [B-scored] LIYXEY_comp_0       13.11s pass=True clash=0 worst=0.7513
  [B-scored] YENDUS_comp_0    KILLED at hard cap 330s
  [B-scored] YENDUS_comp_0           ?s pass=False clash=- worst=-

================ SUMMARY ================
  A-default  pass 16/22  median 13.87s  total 1202.9s  >30s: 8  clash 16 over 1/17 mols (severe 7, worst_overlap min 0.4344 med 0.7588)
  B-scored   pass 18/22  median 5.54s  total 351.3s  >30s: 2  clash 2 over 2/19 mols (severe 0, worst_overlap min 0.7283 med 0.7588)
  PASS REGRESSIONS (A pass -> B fail): none
  PASS FIXES       (A fail -> B pass): ['GAVSED_comp_0', 'QIDKUL_comp_0']
  wrote /tmp/claude-1000/-home-tjmustard-Documents-GitHub-OIN-SMILES/ecbea3d0-9b2a-4420-bb7e-33c571123da0/scratchpad/ab_v2.json
```
