# Round-trip residual tail — R4 triage & fix

Session: `feature/roundtrip-residual-tail`, branched from `main` @ `c06495c`
(R1/R2/R3 landed). Scope: the round-trip long tail no active session owned —
the non-binding donor-H atom-count residual, plus triage/routing of
`string_mismatch_other`, `atom_stereo`, and three singleton buckets.

## Task 0 re-measurement (branch `18198811`, 2026-07-11)

The CASE_REGISTRY the original triage below was written against predated R1/R2/R3
and was regenerated away, so every routed bucket was **re-run and re-classified**
on branch `18198811` (= `origin/main` `c06495c` + the R4 O/S charge fix), rdkit
2026.3.3, `--quick` FF, over the 68 named bucket molecules. Result: **42 success /
26 fail**. Fresh disposition of the 26 residual failures:

| Class (fresh) | n | Molecules | Owner / disposition |
|---|--:|---|---|
| `donor_H_atom_count` | 12 | TIDJIZ; FABPEG FOKDAM WAMWUE XIYJEU ZOFREU (N-analog, +H); ESOSOU HOSXUJ IWAZAJ NOBYOU QOXPAU RUWZIS (atom-**loss**, −H) | out of R4's O/S scope — see below |
| `atom_stereo` | 11 | AJOKUH APACAW_comp_1 CUQVUF FAMFUV GUXPAS ICOLOD IROXAP NONHUU SUNROK XEMSAK YOSYEM | routed → R5/S6-stereo |
| `high_rmsd` | 2 | ACOXOH (1.09), QEGJOE (1.02) | S5-metrics (borderline FF geometry, not E/Z) |
| `geometry_NON` | 1 | DEKQAN | S3 bond perception |

**What this confirms and what changed vs. the pre-R3 triage below:**
- **O/S fix verified: 22/22.** All 21 targeted O-deficit `donor_H` rows **plus
  BEJSUH** (a 138-atom Ni porphyrinoid whose non-binding S-oxide now serializes
  `S([O-])[O-])`) flip failed → success. BEJSUH round-trips byte-identical.
- **`string_mismatch_other` mostly evaporated under R3.** Of the 28, **15 now pass**
  (11 of 13 E/Z rows, CILGEM + EJUBUH of the @-rows, and KAHZEB + KEDLUA). The
  stale E/Z / `H_on_terminal_oxo_imido` / now-passing @-routes were **removed from
  `tools/triage_overrides.json`**; only the 11 still-failing @-rows remain routed.
- **`H_on_terminal_oxo_imido` is empty on these rows.** The 5 previously force-routed
  rows (FENMIX, KENTEE, KUNWIB, TIPCUR, ZIRPEX_comp_1) all now round-trip
  byte-identical (R3's comparator collapses the terminal azide/imine notation) —
  their overrides were removed.
- **TIDJIZ** (the second named donor-H tail row) fails `donor_H_atom_count` 98→**101**
  (+3 H). The heavy-atom skeleton is identical (`s1==s2` except a κ¹-pyrazolyl slot
  permutation); the 3 phantom H land on the **three κ¹-pyrazolyl N donors** of its
  tris(pyrazolyl)borate scorpionate. This is the **aromatic-N analog** of the O/S
  fix — deliberately outside R4's O/S scope (bare/aromatic non-binding N overlaps
  the nitride/ammine notation owned by `oin/inline.py`, R3). Left in
  `donor_H_atom_count`, documented, candidate for the N-notation follow-up.

**No new encoder code was needed:** the landed O/S fix already resolves every O/S
case in these buckets; all remaining failures are out of R4's scope (N-analog
phantom-H, generator atom-loss, embed atom-stereo → R5, borderline FF geometry, and
DEKQAN's bond-perception `g:NON`).

## Fix landed — non-binding 0-H chalcogen donors (`donor_H_atom_count`)

**Mechanism.** `OINSanitizer.generate_robust_smiles` (`utils/oin_aligner.py`)
locked the hydrogen count only for metal *binding* atoms. A **non-binding**
heteroatom with a valence deficit (a croconate/oxo ring O, a nitrito `–O`) has no
H and no metal bond to fill its valence, yet `MolToSmiles` serialized it BARE
(`c(O)`, `ON=O`) because `SetNoImplicit` alone does not force a bracket. The
MetalloGen adapter's `MolFromSmiles` then re-added the implicit H, so the
regenerated 3D structure gained a phantom hydrogen the input never had →
atom-count mismatch (COLWIK croconate 55→58, ACOXEX nitrito 75→77).

**Why charge, not radical.** A neutral 0-H oxygen is a *radical*; the adapter
drops `[O]` back to bare `O` (verified — it stays phantom), and UFF cannot type a
radical (would trade an atom-count fail for `no_conformers`). A **formal charge**
(a valence-1 O → `[O-]`, phenolate/alkoxide) is 0-H, closed-shell, survives the
adapter unchanged, and embeds. The fix charges each non-binding O/S deficit by its
deficit. Encoder-only — **no adapter change needed**.

**Safety.** Restricted to O/S. An aqua/hydroxo/carbonyl/ether O sits at full
valence (no deficit) and a real O–H shows a bonded H, so neither is ever touched.
N is excluded on purpose (it overlaps the nitride/ammine notation ambiguity owned
by `oin/inline.py`).

**Verification.**
- 21/21 targeted O-deficit `donor_H` rows flip **failed → success** (`--quick`,
  FF): ACOWIA, ACOXEX, CIRNEX, COLWIK, EHOSUP, EMEGIM, FEWZIT, IKUMIL, IMIYIO,
  IWOLEN, KAKGAG, KAPTON, NAVPIJ, PAMJOF, PAQQAC, QAJDIR, QAJREB, SOXVAF, VETBEE,
  VIMQEO, VOTYOU.
- Regression: the four S1 ammines (AFAVIO, OQIHUT, RIZVAY, XILBIF) keep NH₃ and
  stay success; a 60-molecule passing sample shows **0** molecules gain `[O-]`
  (blast radius is essentially the failing rows only).
- Full unit suite green on both blessed rdkit versions: 370 OK / 5 skip (2026.3.3),
  375 OK / 5 skip (2025.9.3). New guard `tests/unit/test_non_binding_donor_hydrogens.py`
  (5 tests; 2 fail against pre-fix code).

## `string_mismatch_other` (28) — triaged, all routed (none R4-fixable)

> **Historical (pre-R3 triage).** Superseded by the Task-0 re-measurement above:
> on current code 15 of these 28 pass and only 11 @-rows remain routed. Kept for the
> per-row RDKit-diff reasoning.

Per-fragment RDKit canonicalization (slot/winding/metal-geo tokens stripped) shows
**every** row is stereo-only:

| Class | n | Molecules | Routed to |
|---|--:|---|---|
| E/Z bond direction | 13 | ACOXOH, ACUTAU, BUCFOV, FIXYOB, FOPWUC, HOHTAA, KASYEL, KAZJAY, LUQPET, PITHUY, QEGJOE, VUWTIR, YEJFIC | `EZ_bond_stereo` |
| sp3/@ atom chirality | 13 | AJOKUH, APACAW_comp_1, CILGEM, CUQVUF, EJUBUH, FAMFUV, GUXPAS, ICOLOD, IROXAP, NONHUU, SUNROK, XEMSAK, YOSYEM | `atom_stereo` |
| canonicalize-EQUAL | 2 | KAHZEB, KEDLUA | expected R3-resolved (left unrouted; re-measure should show success) |

Routing entries added to `tools/triage_overrides.json`. (Triage ran on the stored
pre-R3 pairs; R3's comparator likely already collapses several E/Z rows to success —
a full re-measure of these non-owned rows is out of R4 scope.)

## `atom_stereo` (25) — recommend a dedicated **R5-stereo** session

Sampled AHEBEV, BABWAD, FADSAE, YAXVOJ, DAXJUI (plus the string_mismatch stereo
rows above): every one is an sp3-carbon (and some P) stereocenter inversion
(`[C@H]`↔`[C@@H]`), frequently entangled with eta-ring winding markers. The
MetalloGen embed builds the wrong diastereomer/enantiomer because it does not
enforce the input's atom chirality; the encoder round-trips the geometry it is
given. This is `core/chirality.py` (`CIPAssigner`, `ChiralityRecoveryUtility`) +
the generator-embed stereo-enforcement path — a different risk profile from R4's
H-count work, and entangled with the known latent `_apply_double_bond_stereo`
formal-charge bug. **Recommendation:** open an **R5-stereo** session owning
`core/chirality.py` + the embed stereo path (`generator3d/`); it should take the
25 `atom_stereo` rows plus the 13 `string_mismatch_other` @-rows routed above.

## Singletons — routed, not fixed here

- **`geometry_NON` (DEKQAN).** Encoder emitted `[Y_NON]` with three radical
  nitrates `[O]N([O])=O` — the Y–O(nitrate) bonds were never perceived, so the
  coordination number is wrong and no geometry template matched. Root cause is
  bond perception (`utils/xyz2mol.py`, S3), not the encoder/adapter. Route to S3.
- **Atom-count *decrease* rows in `donor_H_atom_count`** (ESOSOU 86→85, HOSXUJ
  33→31, IWAZAJ 44→43, NOBYOU 37→35, QOXPAU 49→47, RUWZIS 63→62). These *lose*
  atoms on generation — a generator/perception atom-loss defect, not encoder
  phantom-H. IWAZAJ/NOBYOU/QOXPAU have `smiles_1 == smiles_2` (INENOF-class,
  downstream in the MetalloGen build → generator/R2); ESOSOU/HOSXUJ/RUWZIS also
  show fragment reordering (comparator/perception). Left in `donor_H_atom_count`
  and documented here rather than force-fit to a wrong class.
- **Non-binding deficit-N rows** (FABPEG, FOKDAM, WAMWUE, XIYJEU, ZOFREU). The
  N analog of the fix above: a non-binding azide/ketenimine/macrocyclic-amine/
  amidine N serializes bare and gains a phantom H. Deliberately outside R4's O/S
  scope — bare non-binding N overlaps the nitride/ammine notation ambiguity owned
  by `oin/inline.py` (R3). Candidate for a follow-up N-notation session.
