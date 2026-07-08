# S6 — @-stereo (P/N donors) + generator E/Z honoring (task #8) + embed corners

Branch: `feature/roundtrip-stereo` · Read `docs/handoffs/README.md` first.

## Mission

Three stereo-flavored residuals, previously deferred:

1. **@-stereo mismatches** (31 rows): input and re-encode differ only in @/@@
   tags — P and N stereocenters whose handedness is lost or flipped through
   generation.
2. **Generator ignores C=C E/Z** (7 rows, "task #8"): the OIN encodes `/C=C\`
   faithfully (encoder fixed in fb8505f), but the 3D generator embeds an
   unconstrained (cis-biased) geometry.
3. **MetalloGen no-conformer corner cases** (triage, 3 fresh rows).

Success = fresh @-stereo and E/Z rows round-trip; no-conformer corners triaged
into fix/document buckets.

## Evidence pack

@-stereo fresh: `QEXKOT_comp_0` [Ni_SPL], `SUXJOL_comp_0` [Rh_SPL],
`WODBIC_comp_0` [Rh_SPL], `WIWYIM_comp_0` [Rh_SPL] — skeletons identical,
only @ tags differ (classifier evidence: "@ tags differ, skeleton identical").
E/Z fresh: `XIZXAG_comp_0` [Zn_TPY]. Known related: `VOacac2` fails in
`tests/integration/verify_xyz_to_oin.py` with spurious `/C=C(/C)` on acac —
same class, use it as a fast local repro.
No-conformer corners: `HIQCIU_comp_0` [Cu_TET] (macrocycle whose m-SMILES
carries `/C=C/` constraints), `KEBBUO_comp_0` [Mo_TET] (spiro-siloxane with
`[Si@@]`/`[Si@]`), `PUMKEN_comp_0` [Re_OCT] (mesitylene σ-CH₃ carbon donor
`[CH3:4]` — can MetalloGen place a CH₃ donor at all?).

## Prior art — READ BEFORE CODING (hard-won conventions)

- **Zone-A chiral P donor fix (commit 65255a1, v0.3.5)**: `build_contract_mol`
  stamps `_OIN_CIPCode_LP` + seeds a `CHI_TETRAHEDRAL_CW` tag for Zone-A P
  donors gated on the parsed `[P@]` tag, AFTER 3D perception, so
  `ChiralityRecoveryUtility.recover()`'s lone-pair verify-and-flip keeps and
  orients it. **CRITICAL: the label must come from `rdCIPLabeler` on the
  metal-free template (`_template_lp_label`), NOT legacy
  `Chem.AssignStereochemistry` `_CIPCode` — the two DISAGREE for 3-coordinate P
  (ACUWUT: legacy 'R' vs rdCIPLabeler 'S') and recover() recomputes with
  rdCIPLabeler.** Sourcing the legacy label round-trips the WRONG enantiomer.
- Backbone P/S/Si stereo (494629c): template `_CIPCode` stamped as
  `_OIN_CIPCode` on backbone P; Si/S in the perceive-then-flip carry set.
  Guards: `tests/unit/test_zone_a_p_donor_stereo.py`,
  `tests/unit/test_backbone_heteroatom_stereo.py`.
- The MetalloGen embed is stereo-blind for lone-pair chirality → handedness is
  random per seed → stereo must be re-asserted from the template, never
  trusted from geometry. Your fresh @ cases are likely paths that never got the
  65255a1 treatment (e.g. N stereocenters, or P in geometries/fragment shapes
  the gate misses) — diagnose which stamp/recover branch drops the tag.
- **Task #8 root cause note (unverified)**: `src/oinsmiles/generator3d/ligand.py:8`
  parses with `sanitize=False` → the internal ace_mol drops bond stereo → embed
  unconstrained. Verify, then thread bond-stereo through to the embed (a
  dihedral constraint or a post-embed flip-and-reoptimize selection, cheapest
  correct thing first).
- Zone-A N is KNOWN-deferred (RDKit clears trivalent `[N@]`; needs an
  out-of-band marker — Option C design). If your @ cases turn out to be Zone-A
  N, document and defer again rather than fighting RDKit.

## Verify-first steps

1. `--only QEXKOT_comp_0` — is the input tag on P or N? Zone-A (metal-bound) or
   backbone? Which branch of `recover()` clears/flips it? (Instrument
   `core/chirality.py` else-branch ~547 and the `_OIN_CIPCode*` stamps.)
2. `--only XIZXAG_comp_0` + `VOacac2` — confirm ligand.py:8 hypothesis by
   checking whether the ace_mol still carries bond stereo after parse.
3. For each no-conformer corner: get the m-SMILES from the error, try
   MetalloGen on it directly, and bin as (a) fixable input massage, (b)
   MetalloGen upstream limitation → document.

## Files

- **Own:** `src/oinsmiles/core/chirality.py`,
  `src/oinsmiles/generator3d/ligand.py`, and if needed the embed call path in
  `generator3d/embed.py` (+ new test files
  `tests/unit/test_donor_stereo_roundtrip.py`,
  `tests/unit/test_generator_ez_honoring.py`).
- **Shared-caution:** if the fix requires touching `build_contract_mol` stereo
  stamps, coordinate with S2 (they own that function this wave) — propose the
  diff in your PR body rather than editing it in parallel.
- **Read-only:** `convert_parsed_to_msmiles` (S1), `utils/xyz2mol.py` (S3),
  `oin_aligner.py`/`compare.py` (S4).
- **Regression floor:** `test_zone_a_p_donor_stereo.py`,
  `test_backbone_heteroatom_stereo.py`, `test_double_bond_stereo_encoding.py`.
  Run stereo tests across ≥3 embed seeds and on BOTH blessed rdkit versions.

## Acceptance

- ≥3 of the 4 fresh @-stereo molecules round-trip with matching tags across
  seeds (any remainder documented with the exact dropped-tag branch).
- `XIZXAG_comp_0` and `VOacac2` round-trip with correct E/Z.
- No-conformer corners: each has a verdict (fixed / documented-unsupported with
  a clear error message).
- Full unit suite green; the flaky task-#8 generator E/Z test noted in the
  worklogs must be made deterministic or its flakiness root-caused as part of
  this work.
