# Documented limitation: carborane / 3c2e borane clusters (no session this wave)

*(v0.4.2 refresh of `spec/handoffs/v0.3.6/wontfix-carboranes.md`. Counts are a dated snapshot of a
live-growing accumulator, so they are quoted as **c7edeeb6-floor / current backlog**, per the wave
protocol: the floor is a set of IDs, never a percentage.)*

The `carborane_unsupported` rows fail in the **FORWARD** encode: `get_lig_mol`
(`src/oinsmiles/utils/xyz2mol.py:426`, which calls `AC2mol` at `:451-462`) cannot build a valid
RDKit template for ligand fragments like `[H]B1[B-]2([H])[B-]3([H])...` — polyhedral boranes and
carboranes whose 3-center-2-electron (3c2e) bonding has no faithful 2-electron SMILES valence model.
`AC2mol` is a **two-centre** bond-order solver; the multi-centre bonding these cages need is outside
the model it implements, and the charge sweep does not rescue it.

The row fails with a specific message, e.g. `HIMQIF_comp_0`:

```
XYZToSMILES failed: ValueError: xyz2mol failed: get_lig_mol failed for ligand
fragment #0 (SMILES: '[H]B1[B-]2([H])[B-]3([H])B(...)[B-]4([H])[B-]([H])...')
```

This is a **notation-design problem, not a bug**: OIN would need an explicit cluster convention
(e.g. treating the cage as an eta-like multi-atom unit or a pseudo-atom) before these can round-trip.
Deferred until a design exists. Interim expectation: the encoder should fail with a clear
"polyhedral borane ligands unsupported" message rather than a raw RDKit valence traceback (S3 may
route the *message*, but the capability itself is out of scope for this wave). Keep the rows under
`wontfix-docs` so they are not mistaken for regressions.

## Members (snapshot 2026-07-14)

**c7edeeb6-floor (36 — reproduce on the pinned baseline commit):**

```
AVOFIB_comp_0 BEKLUA_comp_0 BEKMIP_comp_0 CAKBEW_comp_0 CAKBOG_comp_0 COZCEZ_comp_0
GANYEZ_comp_0 GOHWOQ_comp_0 HAXJAS_comp_0 HAXJOG_comp_0 ICEZIC_comp_0 JABGAX_comp_0
JAFMIP_comp_0 JAFTAO_comp_0 JAFTES_comp_0 MAFSIY_comp_0 MODZUA_comp_0 OZAREO_comp_0
PAQBOZ_comp_0 PAQCAM_comp_0 PAYTUH_comp_0 PEKQII_comp_0 PEKQUU_comp_0 RANCIU_comp_0
RANMUR_comp_0 RAWJEG_comp_0 RIWKAK_comp_0 RIWKEO_comp_0 RONPES_comp_0 RONQET_comp_0
RONQOD_comp_0 RULBUV_comp_0 ULODUU_comp_0 XUKRIF_comp_0 YIBZIV_comp_0 YIVLAQ_comp_0
```

**Current backlog (92 — mixed-provenance, grows as the `--quick --continue` accumulator runs).**
Regenerate any time with `tools/classify_failures.py --output-dir <copy>` on a **copy** of
`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0/` (never the live dir); the class is `carborane_unsupported`.
The additional rows over the floor (as of this snapshot):

```
AFOGEK_comp_0 AQUHOL_comp_0 AQUJIH_comp_0 BAFWAI_comp_0 BAFWOW_comp_0 BAFXAJ_comp_0
BELBEB_comp_0 BIWJID_comp_0 CAJZOD_comp_0 CAKBIA_comp_0 CEMTIZ_comp_0 CEQWAY_comp_0
CEQWEC_comp_0 CEQWOM_comp_0 EGETOB_comp_0 GATHIT_comp_0 GIFVUN_comp_0 GOHWEG_comp_0
GOHWUW_comp_0 GOHXAD_comp_0 GUHMOL_comp_1 GUNZUK_comp_0 HAPCUY_comp_0 HAXJEW_comp_0
HAXJIA_comp_0 HIMQIF_comp_0 HIMQOL_comp_0 ICEZUO_comp_0 JAFSER_comp_0 LETYAL_comp_0
MAFSAQ_comp_0 MAFSEU_comp_0 NEGVOL_comp_0 OLUZIE_comp_0 OTOLAL_comp_0 PAQBUF_comp_0
PEKQOO_comp_0 PEKRAB_comp_0 POVVAZ_comp_0 QAMJAP_comp_0 QIDHIU_comp_0 QODDOF_comp_0
RAJNUO_comp_0 RIWKIS_comp_0 RIWKOY_comp_0 RIWKUE_comp_0 RONQAP_comp_0 RUJSIA_comp_0
SOQKEQ_comp_0 UFUCUU_comp_0 ULOFAC_comp_0 ULOFEG_comp_0 UYESUO_comp_0 UZUYAQ_comp_0
XOJKEM_comp_0 YIVQOJ_comp_0
```
