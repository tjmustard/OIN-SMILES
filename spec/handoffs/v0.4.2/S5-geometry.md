# ▶ START HERE — S5 geometry (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** S5 is the
**merged geometry + winding phase** — `winding_flip` folds in here because its winding-aware
conformer pick lives *inside* `_select_by_geometry` and is data-coupled to that function's `scored`
list (verified `metallogen_adapter.py:1180`). One owner avoids a mid-function collision.

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-geometry -b feature/roundtrip-geometry release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-geometry && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- prior — `spec/handoffs/v0.3.6/S4-eta-winding.md`, the CN-8 `SQA` precedent (`tests/unit/test_cn8_geometry.py`)
- `tools/triage_overrides.json` — its `_comment_winding` explains why a differing winding multiset is
  a real coordinated-face difference (encoder canonicalizes reversible-ring winding), and lists rows
  (QOFTOU) that are generator stereochemistry, not notation.

### 3 · Verified code paths
- **geometry_or_fragment_change (53):** ENCODER borderline classifier + GENERATOR energy fallback.
  - Classifier `oin_aligner.py:1489 classify_coordination_geometry` → `:796 _find_best_geometry_match`
    → `:805 _match_geometry_candidates`. For CN4 it scores **SPL/TET/TPY** and keeps strict min-RMSD
    (`:865`, **no tolerance band** — a near-tie flips the label). `:874 _map_to_template` fits
    **unit vectors only** (purely angular), so a flattened tetrahedron scores closer to TPY than TET.
  - Generator `metallogen_adapter.py:1121 _select_by_geometry` filters conformers re-classified as the
    target (`:1164-1166`) but returns `mols[0]` (lowest energy, any geometry) when none match
    (`:1199-1213`).
- **geometry_NON (4, e.g. DEKQAN [Y_NON]):** ENCODER — no template for **CN≥9**. `oin_aligner.py:700`
  returns `g:NON`; candidate table falls back to `["OCT"]` for `n>8` (`:840-841`) → `n>len(vectors)`
  `continue` (`:858-859`) → `best_result=None`. Adapter rejects at `metallogen_adapter.py:114-116`;
  geo dict has no CN-9 key (`:73-89`). (CN-8 `SQA` already landed in v0.4.1.)
- **winding_flip (29, e.g. AGOVOK):** GENERATOR geometry — the encoder is already correct (fixed `>`
  for reversible rings `oin_aligner.py:1375-1404`; geometric sign only for load-bearing rings via
  `:1410 _determine_winding` / `oin/winding.py:signed_circulation`). The generator's winding-aware
  pick `metallogen_adapter.py:1173-1197` (consumes `scored` at `:1180`) failed to sample a conformer
  matching the target multiset and fell back to the wrong face; pool widening `:1264-1275`
  (`ETA_SELECT_POOL`), helper `_eta_winding_multiset :1043`.

### 4 · Mission & scope guard
- **geometry classifier:** add a **TET/TPY hysteresis margin** so a near-tie prefers the input's
  label / does not flip on a small angular perturbation (`_match_geometry_candidates :865`). Prove it
  does not regress genuine TPY/SPL cases.
- **geometry_NON:** add a **CN≥9 template** to `TEMPLATE_SPECS` + the candidate table (mirror the
  CN-8 `SQA` landing) and a matching key in the `OIN_TO_METALLOGEN_GEO` **dict** (`:73-89`). Then the
  `:114-116` raise stops firing **without editing it** — leave `:114-116` alone (S1's function).
- **winding_flip:** strengthen `_select_by_geometry`'s winding pick — widen/vary the conformer pool
  (`ETA_SELECT_POOL`, `:1264-1275`) so a load-bearing allyl/diene face is actually sampled, or add a
  geometric winding *construction* rather than relying on re-encode sampling. **Re-triage first**:
  some winding rows are misrouted generator stereochemistry (QOFTOU builds rac/meso
  non-deterministically) — those may be irreducible; document them.
- Because winding shares `_select_by_geometry` with geometry selection, keep both changes coherent —
  the `scored` list you feed the winding loop must remain geometry-first.

### 5 · Owned files (edit only these regions)
- `src/oinsmiles/utils/oin_aligner.py` — geometry matcher/classifier/templates
  (`_match_geometry_candidates :805`, `_map_to_template :874`, `TEMPLATE_SPECS`,
  `classify_coordination_geometry :1489`; add CN≥9 + hysteresis). The winding **canonicalization**
  (`:218,:1375-1410`) is correct — read-only; don't change encoder winding.
- `src/oinsmiles/generation/metallogen_adapter.py` — **`_select_by_geometry :1121-1213`** (incl.
  winding block `:1173-1197`) + **`generate` pool `:1264-1275`** + **`_eta_winding_multiset :1043`**
  + **`OIN_TO_METALLOGEN_GEO` dict `:73-89`** (dict only). Do **NOT** touch `convert_parsed_to_msmiles`
  `:105-238` (S1) or the sp3 stamp `:559,:776-778` (S6b).

### 6 · Gate
- geometry-change goldens: input and re-encode agree on the geometry label; NON goldens generate
  (no raise) and round-trip; winding goldens: the re-encoded multiset matches the target.
- **Non-`--quick`, `--mol-timeout 1800`** for winding/geometry goldens (conformer-yield sensitive).
- **No CN regressions** — per-CN spot-check (a CN4 fix must not regress CN5/CN6; run the four v0.4.0
  goldens cisplatin/ferrocene/Ir(ppy)3/BINAP and confirm byte-identical or better).
- New guard tests under `tests/unit/` (CN≥9 template; TET/TPY hysteresis; winding sampling), failing
  pre-fix. Full unit suite green (own baseline first); `ruff check` clean.

### 7 · Landing
Squash-merge into `release/v0.4.2` (see `SESSION_PROMPTS.md`).
