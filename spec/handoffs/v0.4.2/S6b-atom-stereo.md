# ▶ START HERE — S6b atom stereo + [S@SP3] (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** S6b is
function-disjoint from S6a in `embed.py` (you own `_apply_atom_chirality`/`_permutation_is_odd`; S6a
owns `_apply_double_bond_stereo`) and owns `core/chirality.py` alone — run concurrently. **Only if
you need the `chiral_centers` capture (`ligand.py:73-83`) do you chain after S1** (audit that first).

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-atom-stereo -b feature/roundtrip-atom-stereo release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-atom-stereo && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- prior — `spec/handoffs/v0.3.7/R5-atom-stereo.md`, `docs/KNOWN_LIMITATIONS.md` ("sp3/heteroatom
  atom-stereo" — R5 took this 0/25 → 18/25; the residuals are documented)
- **`tools/triage_overrides.json`** — RE-TRIAGE YOUR INPUTS HERE FIRST (see §4)

### 3 · Verified code paths
- **atom_stereo (18, e.g. CUQVUF):** residual after R5, BOTH sides. Encoder
  `core/chirality.py:367 assign_all` (tags from 3D `:403`), `:504 recover` (sp3 re-orient branch
  `:562-598`, keyed on `_SP3_CIP_PROP`), `:50 _reparse_aromatic_cip_label` (metal-free fresh
  re-parse, kekulize-skipped `:63-73`). Generator `embed.py:480 _apply_atom_chirality`,
  `:470 _permutation_is_odd`; metal-free template label stamped
  `metallogen_adapter.py:559 _template_sp3_label` / `:776-778`. CUQVUF-type inversions are where the
  contract-mol flip loop and the fragment CIP disagree and the `_SP3_CIP_PROP` re-orientation misses
  the centre (e.g. P/S not in `(6,14,16)` at `chirality.py:575`, or no legacy `_CIPCode`).
- **`[S@SP3]` string_mismatch subset (e.g. BAZMEX `CCS{4}` vs `CC[S@SP3]{4}`):** ENCODER — the
  sulfur donor's geometry-derived tetrahedral tag is never cleared. Origin
  `chirality.py:403 AssignAtomChiralTagsFromStructure`; the clear at `:307
  _clear_spurious_high_coordination_stereo` handles only OCT/TBP (`:325-328`), and `recover()`'s
  degree-keyed cleanup (`:604-605`) bails on anything not in `_PN_ATOMIC_NUMS={7,15}` — **S(16) is
  never reached**, so its tag survives into the SMILES (emitted at `oin_aligner.py:356` /
  `xyz2mol.py:1219`).

### 4 · Mission & scope guard
- **RE-TRIAGE FIRST.** v0.3.6 found **12 rows misfiled as `atom_stereo`** that were η-hapticity slips,
  geometry-class changes, or a destroyed SF₅. And some `string_mismatch_other` rows are `@`-stereo
  while others are η-slips. Reclassify via `tools/triage_overrides.json` (never by hand): misfiled
  geometry/winding rows route to **S5**, not S6b. Verify stereo **only through
  `XYZToSMILES.convert()`** — `get_tmc_mol` defaults `with_stereo=False` and fabricates false passes.
- **`[S@SP3]` (encoder):** extend the clear to sulfur — either add S(16) handling to
  `_clear_spurious_high_coordination_stereo :307` (strip a `CHI_TETRAHEDRAL` tag on a donor S with no
  supporting `_OIN_CIPCode*`), or widen `recover()`'s degree-keyed branch (`:604-634`) beyond
  `_PN_ATOMIC_NUMS`. Do not strip a genuine diastereomer — gate on "no supporting CIP property".
- **atom_stereo (both sides):** extend the `_SP3_CIP_PROP` re-orientation to the centres R5 missed
  (widen the atomic-number gate at `chirality.py:575` where the chemistry warrants; ensure the
  metal-free fresh re-parse `_reparse_aromatic_cip_label` is used on **both** the stamp and the
  comparison — using it on one side flips a P-and-arene-adjacent centre). Route the genuine
  "resolved-stereo-but-fails-another-class" rows (donor-H, high_rmsd, no_conformers) to their owning
  phase — they leave `atom_stereo`.

### 5 · Owned files (edit only these regions)
- `src/oinsmiles/core/chirality.py` — the clear/recover/reparse machinery.
- `src/oinsmiles/generator3d/embed.py` — **`_apply_atom_chirality :480` / `_permutation_is_odd :470`
  only**. Do **NOT** touch `_apply_double_bond_stereo :417` (S6a) or `get_embedding :598` (S7).
- `src/oinsmiles/generation/metallogen_adapter.py` — **the sp3 CIP stamp `_template_sp3_label :559`
  / `build_contract_mol` `:776-778` only**. Do **NOT** touch `convert_parsed_to_msmiles :164-222`
  (S1), `_select_by_geometry :1121-1213` (S5), or the geo dict `:73-89` (S5).
- **Audit-only** `ligand.py:73-83` (`chiral_centers` capture): if you must change it, chain after S1
  and coordinate — it is index-coupled to S1's `AddHs`.

### 6 · Gate
- atom_stereo goldens (post-re-triage) + `[S@SP3]` goldens round-trip: `@`/`@@` and donor tags match,
  canonical key matches, via `--only` into `/tmp/rt-atom-stereo`, measured through
  `XYZToSMILES.convert()`.
- New guard tests under `tests/unit/` (S-tag clear; the missed sp3 re-orientation), failing pre-fix.
  Extend `tests/unit/test_heteroatom_atom_chirality.py`.
- Full unit suite green (own baseline first); `ruff check` clean. Currently-passing R5 cases
  (KAPCEM-class, the 18 successes) must **stay** passing.

### 7 · Landing
Squash-merge into `release/v0.4.2` (see `SESSION_PROMPTS.md`).
