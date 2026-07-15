# v0.4.2 BASELINE — clean single-commit floor

Baseline commit **`c7edeeb6`** (= tag `v0.4.1`). Built by P0 **without a new sweep**: the live accumulator had already produced a single-commit sample on this exact commit, so the floor is extracted directly (provenance filter + confirmed-on-baseline goldens). It grows as the accumulator runs — regenerate with `scratchpad/build_baseline.py`.

## Must-not-regress passing set (trusted provenance)

**5964 molecules** stamped `commit_id == c7edeeb6` and `status == success`. Full ID list: `spec/handoffs/v0.4.2/baseline_pass_c7edeeb6.txt`. The capstone's blocker gate is `{{passes on release/v0.4.2}} ⊇ {{this set}}`.

**Untrusted passers (9662):** stamped on the dirty/older tree (`5538b722-dirty` etc.), not `c7edeeb6` — a pass there may not hold on the baseline. Treat as *unverified*; do not add to the must-not-regress set without a baseline re-run.

## Per-class floor (failures confirmed on `c7edeeb6`)

Counts are failures **already stamped `c7edeeb6`** (the clean floor sample so far); the wider mixed-provenance backlog is larger. Goldens are drawn from these baseline-stamped rows, so each reproduces on the exact baseline commit. Full per-class member lists: `spec/handoffs/v0.4.2/baseline_fail_c7edeeb6.json`.

### Fixable classes (wave targets)

#### `donor_H_atom_count` — 82 on baseline
- goldens: `AJODEI_comp_0`, `APAGEE_comp_0`, `ARONEA_comp_0`, `ATECUZ_comp_0`, `BAXWAA_comp_1`, `BIBROW_comp_0`, `BOLDUC_comp_0`, `BOXJUU_comp_0`
- repr `AJODEI_comp_0`: error: Atom count mismatch at FF_reroll_5. Input 97 != Gen 95

#### `H_on_terminal_oxo_imido` — 2 on baseline
- goldens: `UBUMUC_comp_0`, `YOVMUR_comp_0`
- repr `UBUMUC_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Co_OCT].CC(N{0}=O)=C(C)N{1}=O.CC(=N{2}O)C(C)=N{3}O.CC(C)(C)c1cc(/C=N/c2cccn{4}c2)c(O)c(C(C)(C)C)c1.[CH2]{5}CCC#N, Got: [Co

#### `geometry_NON` — 1 on baseline
- goldens: `XERTUK_comp_3`
- repr `XERTUK_comp_3`: error: Generation/Verification failed at UFF_1: ValueError: Geometry code 'NON' not supported by MetalloGen mapping.

#### `geometry_or_fragment_change` — 29 on baseline
- goldens: `ALAKII_comp_0`, `AXIYOX_comp_2`, `BACCEP_comp_0`, `BIJYAV_comp_0`, `CEBKIE_comp_0`, `DADXOW_comp_0`, `DAJPAH_comp_0`, `DANXEY_comp_0`
- repr `ALAKII_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Ni_SPL].COc1ccccc1[C@H](CCO)[C@@H](N{0}=C(c1ccccc1)c1ccccc1N{1}C(=O)CN{2}1CCCCC1)C(O{3})=O, Got: [Ni_TET].COc1ccccc1[C@H](

#### `winding_flip` — 14 on baseline
- goldens: `ATUROX_comp_0`, `BOXKIJ_comp_0`, `FEQLUI_comp_0`, `FICCIE_comp_0`, `HEHHOU_comp_0`, `IRAYAB_comp_0`, `LOFMAV_comp_0`, `MOQKIN_comp_0`
- repr `ATUROX_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Zr_TET].Cc1c(C)c(C)c{0}2c{0}(CCc{1}3c{1<}(C)c{1}(C)c{1}4c(C)c(C)c(C)c(C)c{1}34)c{0<}(C)c{0}(C)c{0}2c1C.[Cl]{2}.[Cl]{3}, Go

#### `EZ_bond_stereo` — 23 on baseline
- goldens: `AHAZOZ_comp_0`, `AROHIA_comp_0`, `AYUYIE_comp_0`, `BOCGEH_comp_0`, `DERLEU_comp_0`, `DOXRUD_comp_0`, `FOJHES_comp_0`, `GIQWAG_comp_0`
- repr `AHAZOZ_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Cr_TET].Cc{0}1[cH]{0}[cH]{0>}[cH]{0}c{0}(/C=[N+](\[O-])Cc2ccccc2)[cH]{0}1.C{1}#O.C{2}#O.C{3}#O, Got: [Cr_TET].Cc{0}1[cH]{0

#### `atom_stereo` — 11 on baseline
- goldens: `JEKQAS_comp_0`, `JUCCUH_comp_0`, `KEBBUO_comp_0`, `ORIHUU_comp_0`, `POYJIX_comp_0`, `REPZUJ_comp_0`, `SEMTOV_comp_0`, `VEJXOZ_comp_0`
- repr `JEKQAS_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Pd_SPL].Cc1ccc([S@@](O{0})(=O)=N{1}c2ccc3ccccc3c2-c2c(N3C{2}N(C)C=C3)ccc3ccccc23)cc1.[I]{3}, Got: [Pd_SPL].Cc1ccc([S@](O{0

#### `encode_crash_other` — 4 on baseline
- goldens: `ASISAX_comp_0`, `IROXET_comp_0`, `SUNXAB_comp_0`, `XEVMAN_comp_0`
- repr `ASISAX_comp_0`: error: XYZToSMILES failed: ValueError: xyz2mol failed: get_lig_mol failed for ligand fragment #0 (SMILES: '[H]C1C([H])C([H])C2C(C1[H])C1C3NC(C([H])C([H])C4NC5C6NC(C([H

#### `kekulize_encode_crash` — 7 on baseline
- goldens: `JOTJEK_comp_0`, `KAXVOX_comp_0`, `KAXWAK_comp_0`, `LEZWAO_comp_0`, `TIYWUV_comp_0`, `ZENZAW_comp_0`, `ZENZAW_comp_1`
- repr `JOTJEK_comp_0`: error: XYZToSMILES failed: ValueError: xyz2mol failed: cannot kekulize molecule and found no quinoid ring to relax: Can't kekulize mol. Unkekulized atoms: 2 3 7 8 9 11

#### `macrocycle_perception` — 9 on baseline
- goldens: `FIMXIK_comp_0`, `HUTCIK_comp_0`, `JOBYIK_comp_0`, `KOSQOA_comp_0`, `LOLROW_comp_0`, `WARQEO_comp_0`, `WEVNUJ_comp_0`, `XIZXAG_comp_0`
- repr `FIMXIK_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Fe_OCT].CC1=N{0}O[B@]2(c3ccncc3)ON{1}=C(C)C(C)=N{2}O[B@](c3ccncc3)(ON{3}=C1C)ON{4}=C(C)C(C)=N{5}O2, Got: [Fe_OCT].CC1=N{0}

#### `garbled_aromatic` — 2 on baseline
- goldens: `DIXXIS_comp_0`, `ROJQOY_comp_0`
- repr `DIXXIS_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Zn_SPY].c1ccn{0}cc1.Ic1ccc(C2=C3C=CC(=N{1}3)C(c3ccc(I)cc3)=C3C=CC(=N{2}3)C(c3ccc(I)cc3)=C3C=CC(=C(c4ccncc4)C4[CH][CH]C2N{3

#### `string_mismatch_other` → `[S@SP3]` subset — 5 on baseline (S6b target)
- goldens: `BAZMOH_comp_0`, `CIDDAU_comp_0`, `HUGSEI_comp_0`, `LUSKIV_comp_0`, `YUMPIH_comp_0`
- repr `BAZMOH_comp_0`: error: String mismatch at FF_reroll_5. Exp: [Ru_OCT].[Cl]{0}.CCS{1}CCN{2}(C)Cc1ccccn{3}1.c1ccc(P{4}(c2ccccc2)c2ccccc2)cc1.[Cl]{5}, Got: [Ru_OCT].[Cl]{0}.CC[S@SP3]{1}CC

### Artifact classes — context only (NOT part of the floor; S7/docs own them)

- `timeout`: 339 on baseline
- `high_rmsd`: 36 on baseline
- `carborane_unsupported`: 36 on baseline
- `no_conformers`: 115 on baseline
- `gen_exception_other`: 24 on baseline
- `rmsd_mapping_failed`: 1 on baseline

### Other failed classes on baseline (triage)
- `string_mismatch_other`: 20

## Method / caveats

- **No headline percentage.** The floor is a *set of molecule IDs*, per the wave protocol.
- Classes were derived with `tools/classify_failures.py::classify` + `tools/triage_overrides.json` (same routing the registry uses), over `results-v0.4.0/individual_reports/*.json` filtered to `commit_id == c7edeeb6`.
- Artifact classes (timeout/high_rmsd/carborane/no_conformers) are context only — highest cost, lowest diagnostic value; S7 triages them at full fidelity, docs records the residual.
