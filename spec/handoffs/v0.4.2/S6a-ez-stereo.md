# ▶ START HERE — S6a E/Z bond stereo (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** S6a **rebases onto
S1** — you both edit `generator3d/ligand.py::get_ligand_from_smiles`. S1 owns the `AddHs` line
(`:63`); you own the `near_donor` block (`:55-62`) immediately above it. **Confirm S1 has landed on
`release/v0.4.2` and rebase before you start**, so your index frame matches.

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-ez-stereo -b feature/roundtrip-ez-stereo release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-ez-stereo && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- the guard that documents the latent bug — `tests/unit/test_chelate_locked_ez.py`
  (`test_dangling_imine_is_not_enforced_pending_charge_fix`); the warning comment
  `generator3d/ligand.py:44-54`
- prior — `spec/handoffs/v0.3.7/R1-charge-aware-contract.md`, the R4 residual triage in
  `docs/roundtrip_residual_triage_R4.md`

### 3 · Verified code paths (this is the documented latent bug)
- **The bug** — `src/oinsmiles/generator3d/embed.py:417 _apply_double_bond_stereo`; the force-set at
  `:450-457` promotes a carried bond to `DOUBLE` and only backs out on a **valence** failure — it
  **never adjusts formal charges**. On AFECIZ the PuLP charge assignment wants that C=N *single* with
  a charged N, so `SanitizeMol` rejects a 4-valent N and **every `ff_clean` raises** (553 s → full
  250-attempt budget at ~1.6 s each).
- **The intentional mismatch that keeps it latent** — `generator3d/ligand.py:55-62` (`near_donor`
  proxy) is deliberately *broader* than the encoder's chelate-ring test
  (`core/translator.py:12 _clear_chelate_locked_bond_stereo`, called `:112`), so the encoder emits an
  E/Z the generator does not enforce (a free C=N on a monodentate arm: AFECIZ, XIZXAG, and the
  N-oxide `/C=[N+](\[O-])` class). The comment at `ligand.py:44-54` says: narrowing `near_donor` to
  the ring test "**detonates a latent bug in `_apply_double_bond_stereo` ... Fix that first, then
  narrow.**"

### 4 · Mission
1. **Make `_apply_double_bond_stereo` charge-aware.** When promoting a bond to `DOUBLE` would
   over-fill a neighbour's valence, adjust the formal charge consistently with the PuLP assignment
   (or decline the promotion and keep the carried single-bond-with-charge), instead of blindly
   force-setting DOUBLE and only checking valence. The bond's E/Z must survive without producing a
   4-valent neutral N.
2. **Then narrow `near_donor` (`ligand.py:55-62`)** from the broad donor-neighbour proxy to the
   encoder's chelate-ring predicate (the one `core/translator._clear_chelate_locked_bond_stereo`
   uses), so a free monodentate C=N (AFECIZ/XIZXAG/N-oxide arms) is **enforced** by the generator and
   round-trips. Only do this **after** step 1 — verify the AFECIZ `ff_clean` no longer raises.
3. `core/translator.py::_clear_chelate_locked_bond_stereo` is **correct** — verify, don't change it.

### 5 · Owned files (edit only these regions)
- `src/oinsmiles/generator3d/embed.py` — **`_apply_double_bond_stereo :417-469` only**. Do **NOT**
  touch `_apply_atom_chirality :480`/`_permutation_is_odd :470` (S6b) or `get_embedding :598` and the
  stereo call-sites `:669,:764` (S7).
- `src/oinsmiles/generator3d/ligand.py` — **the `near_donor` block `:55-62` only**. Do **NOT** edit
  the `AddHs :63` (S1) or `chiral_centers :73-83` (S6b).
- `src/oinsmiles/core/translator.py::_clear_chelate_locked_bond_stereo` — verify only.

### 6 · Gate
- E/Z goldens (AHAZOZ N-oxide, AFECIZ, XIZXAG, the free-arm C=N set) round-trip: the `/`\``
  direction matches, canonical key matches, via `--only` into `/tmp/rt-ez-stereo`. Run `smiles_1`
  diffs under **one** pinned rdkit (the two blessed versions disagree on `/`\`` direction).
- AFECIZ `ff_clean` no longer raises — confirm it embeds within budget (was exhausting 250 attempts).
- Update `tests/unit/test_chelate_locked_ez.py`: the "pending charge fix" xfail/guard now **passes**;
  add a charge-aware `_apply_double_bond_stereo` unit test failing pre-fix.
- Full unit suite green (own baseline first); `ruff check` clean. Currently-passing chelate-locked
  E/Z (VOacac2 / porphyrinoid meso bridges) must **stay** passing — the narrowing must not
  re-introduce spurious markers on ring-locked bonds.

### 7 · Landing
Squash-merge into `release/v0.4.2` (see `SESSION_PROMPTS.md`).
