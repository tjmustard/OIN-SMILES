# S1 — Bare anionic donor protonation + unbound-fragment IndexError

Branch: `feature/roundtrip-donor-h` · Read `docs/handoffs/README.md` first.

## Mission

Generated 3D structures gain hydrogens that the input never had: bare X-type
donors (amide/silylamide/azide N⁻, alkoxide O⁻, acetylide C⁻) come back
protonated, failing the round trip on **atom count** (48 rows) or on
`=O`→`=[OH]` **string mismatch** (28 rows). Separately, the same function
crashes with **IndexError** for fragments that have no binding slot (11 rows).
Success = the evidence molecules below round-trip; both defect classes reach ~0
in the post-fix registry.

## Root cause — PROVEN for the N cases (code read, verify the rest)

`convert_parsed_to_msmiles` in `src/oinsmiles/generation/metallogen_adapter.py`
(donor-H reconciliation block, ~lines 155–200): a bare (non-bracket) binding
atom is only stripped of implicit H when:

- `C`: has a triple bond, is σ-aryl (non-haptic), or `heavy >= 2`
- `O`/`S`: always stripped
- `N`: **only when `heavy >= 2`** ← the bug

A terminal anionic N donor has 0–1 heavy neighbours (bare `N{n}` amide/nitride,
`C[Si](C)(C)N{n}` silylamide, `N{n}N#N` azide) → falls through → RDKit fills
implicit valence → MetalloGen builds NH₂/NH₃.

Key convention to verify then rely on: **the forward encoder brackets every
neutral N–H donor it keeps** (`[NH2]`, `[NH]` — see Track-B rule in the same
block's comment). If that holds (check `OINSanitizer.generate_robust_smiles` /
force-bracket logic in `utils/oin_aligner.py` — read-only for you), then a bare
`N{n}` in an OIN ALWAYS means 0-H amido/imido and the `heavy >= 2` gate should
simply go away for N donors.

## Evidence pack (all reproduced on v0.3.5, rdkit 2025.09.3)

| molecule | defect | detail |
|---|---|---|
| `WAYHOW_comp_0` [Cr_TET] | atoms 62→65 | `N{3}` bare nitride/amide → +3H. smiles_1 == smiles_2 (!) so only the atom-count guard catches it |
| `UDIVUY_comp_0` [Mo_TBP] | atoms 129→131 | `N{4}c1c(Cl)cccc1Cl` anilide → re-encodes `[NH2]{4}` |
| `XADYAC_comp_0` [V_TBP] | atoms 37→39 | `C[Si](C)(C)N{1}` silylamide → `[NH2]{1}` |
| `FENMIX_comp_0` [Mn_OCT] | string | azide `N{0}N#N` → gen `NN#[NH2]{0}` (protonated + wrong-end perception follows) |
| `KIZQER_comp_0` [V_SPY] | string | `N(O)=O` nitro → gen `N(O)=[OH]` |
| `COLWIK_comp_0` [Zn_OCT] | atoms 55→58 | polyphenolate O donors → +3H (check why the O-strip didn't apply — deprotonated NON-binding O?) |
| `INENOF_comp_0` [Au_LIN] | atoms 58→60 | acetylide `C{1}#Cc1ccccc1` → +2H **despite** the triple-bond strip — a second leak, likely elsewhere (embed/ace_mol H fill); verify before assuming this block |
| `TAJBOY_comp_0` [Ti_TBP] | atoms 117→118 | +1H, triage |
| `XUWHOO` (stale row) | atoms 67→66 | nitroso, known Track-B class |

Note the pattern in COLWIK/WAYHOW: the H can land on a NON-binding atom of the
ligand (e.g. the other phenolate O), so fixing only the binding-atom rule may
not clear every case — the OIN template's H counts for the whole fragment are
authoritative; the generated mol must not exceed them.

### IndexError (same function, still live — 11 rows through Jul-8 morning)

```
File ".../generation/metallogen_adapter.py", line 198, in convert_parsed_to_msmiles
  np.argmin(np.linalg.norm(metallogen_vectors - np.array(frag_vectors[0].vector), axis=1))
IndexError: list index out of range
```

`frag_vectors` is empty for a fragment with no binding-slot vector (outer-sphere
counterion / uncoordinated solvent that the encoder still emitted as a
fragment). Cases: `XAXZIH_comp_0`, `PEGXOP_comp_0`, `SOLYEZ_comp_0`,
`CUBDOT_comp_0`, `WAHXOV_comp_1`, `TIGDAO_comp_0`, `NECCIH_comp_0`,
`OLAYUV_comp_0`, `NASZOY_comp_0`, `YUMBEP_comp_0`, `XIFSIR_comp_0`.
Decide the policy explicitly: skip the fragment with a warning (structure then
fails atom-count honestly) vs. raise a clear "uncoordinated fragment
unsupported" error. Look at what the comparator does with uncoordinated
fragments (`normalize_oin_for_comparison` drops empty frags) before choosing.

Also assigned here: 2 stale `geometry_NON` rows (encoder emitted `g:NON`, no
CN template matched — CN-8 SQA landed in v0.3.5; check what CN these are and
whether a template or a clean error is right).

## Verify-first steps

1. Repro: `--only WAYHOW_comp_0,UDIVUY_comp_0,XADYAC_comp_0,FENMIX_comp_0,KIZQER_comp_0`
   — confirm the +H counts.
2. Print the m-SMILES the adapter builds (it appears in "failed to generate"
   errors; add a debug log locally) — confirm the donor carries H before embed.
3. Confirm the encoder force-bracket convention for neutral N–H donors, then
   change the strip rule; re-run step 1.
4. For INENOF (acetylide), instrument where the extra H enters — do NOT patch
   this block blind.

## Files

- **Own:** `src/oinsmiles/generation/metallogen_adapter.py` — ONLY
  `convert_parsed_to_msmiles` (and small helpers you add near it); new test file
  `tests/unit/test_bare_donor_hydrogens.py` (extend the existing
  `test_msmiles_donor_hydrogens.py` patterns; don't rename it).
- **Read-only:** the rest of `metallogen_adapter.py` (S2 owns the template
  functions), `utils/oin_aligner.py` (S4), `utils/xyz2mol.py` (S3).

## Acceptance

- Evidence molecules above round-trip (atom counts equal, keys match) — except
  any you prove to belong to another class (document + reassign in PR body).
- Guard tests: bare N (0/1/2-heavy), azide, silylamide, nitro-O, acetylide,
  plus the dative counter-cases that must KEEP H (`[NH2]`, `[NH]`, `[OH2]`).
- Regression: `ACAWOR_comp_0` and `ABESAD_comp_0` (Track-B carbene/amine fixes)
  still pass; full unit suite green.
