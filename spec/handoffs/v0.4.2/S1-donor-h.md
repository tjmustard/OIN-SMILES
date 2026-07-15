# ▶ START HERE — S1 donor-h (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** S1 **lands first
among the `ligand.py` co-editors** (S6a rebases onto it) because it owns the `AddHs` line that S6a
and S6b's index frames depend on.

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-donor-h -b feature/roundtrip-donor-h release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-donor-h && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- prior diagnosis — `docs/KNOWN_LIMITATIONS.md` (donor-H notation limit), `spec/handoffs/v0.3.6/S1-donor-h.md`
- **do not touch** `_select_by_geometry`, `oin_aligner.py`, `chirality.py`, `xyz2mol.py` (other phases)

### 3 · Verified code paths
- **Encoder bracket convention** — `src/oinsmiles/oin/inline.py:265-317` (`replace_map`): `NH3`→bare
  `N{n}` (`:285`); 0-H nitride gated on `GetDegree()==0`→`[N]{n}` (`:300-306`); pure-organic
  de-bracketed (`:312-315`).
- **Generator H reconcile** — `src/oinsmiles/generation/metallogen_adapter.py:164-222`
  (`convert_parsed_to_msmiles`): only **bare** atoms reinterpreted (`:179`); N branch `strip = heavy
  >= 1` (`:207-219`); chalcogen strip `:203-206`. The nitride/ammine ambiguity is documented at
  `:214-219` ("**ammine reading wins**").
- **H materialization** — `src/oinsmiles/generator3d/ligand.py:63` `Chem.AddHs(..., explicitOnly=False)`.
- **Harness compare** — `tools/test_dataset_roundtrip.py:177-184` ("Atom count mismatch at ...").
  Classifier: `tools/classify_failures.py:131-133` (donor_H), `:153-154` (oxo/imido, a **string**
  mismatch `=\[[ON]H\d?\]` present in gen not input).

### 4 · Mission & scope guard

Two classes: `donor_H_atom_count` (168) and `H_on_terminal_oxo_imido` (9).

**First, sub-distribute the 168** (this class is partly a *notation limit*, not one bug):
- **Fixable (generator strip):** rows where a donor whose H-count the encoder pinned still gets
  wrong `AddHs` — a missing/incorrect `strip` rule in `metallogen_adapter.py:183-219`. Fix these.
- **Notation-limited:** the bare `heavy==0` N where nitride `[N]` and ammine `[NH3]` both serialize
  to `N{n}` and "ammine reading wins" (+3, e.g. AJIJUY 57→60). Genuinely ambiguous without an OIN
  format change. **Route these to `docs`** unless you add a disambiguating encoder marker in
  `inline.py` (extend the `GetDegree()==0` nitride marker at `:300-306` to more cases) — do that
  only if it round-trips without regressing the ammine cases. Promise the subset, not 168.

**`H_on_terminal_oxo_imido` (9, e.g. BADKAS):** generator over-protonation — a **non-binding**
terminal oxo/imido (`N{5}O`, `N{1}=O`) gets an implicit H at `ligand.py:63` because the strip block
only touches the **binding** atom. Add a strip/`SetNoImplicit` rule for terminal oxo/imido/nitroso
neighbours in `metallogen_adapter.py:164-222` (or a valence guard before `ligand.py:63`). The
encoder side is correct.

### 5 · Owned files (edit only these regions)
- `src/oinsmiles/oin/inline.py` (`replace_map` `:265-317`).
- `src/oinsmiles/generation/metallogen_adapter.py` — **`convert_parsed_to_msmiles` only** (`:105-238`,
  H-strip `:164-222`; the geometry `:114-116` raise is yours too but leave it unless geometry_NON
  work overlaps — coordinate with S5, which owns the geo **dict** `:73-89`).
- `src/oinsmiles/generator3d/ligand.py` — **the `AddHs` line `:63`** and any strip you add around it;
  do **not** edit the `near_donor` block `:55-62` (S6a) or `chiral_centers` `:73-83` (S6b).

### 6 · Gate
- Your donor-H fixable goldens + all 9 oxo/imido goldens round-trip (atom count matches; canonical
  key matches) via `--only` into `/tmp/rt-donor-h`, on `c7edeeb6`-based worktree.
- Notation-limited residual explicitly listed for `docs`.
- New guard tests under `tests/unit/` (e.g. `test_terminal_oxo_imido_no_H.py`), failing pre-fix.
- Full unit suite green (measure your own baseline first); `ruff check` clean.
- Currently-passing spot-check: named ammine/amido passers still round-trip (no over-strip).

### 7 · Landing
Squash-merge into `release/v0.4.2` (see `SESSION_PROMPTS.md`). **Announce** so S6a rebases onto your
`AddHs` frame.
