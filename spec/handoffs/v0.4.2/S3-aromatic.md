# ▶ START HERE — S3 aromatic-perception (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** S3 is **fully
file-disjoint** (`utils/xyz2mol.py` + `utils/aromaticity.py`) — run it concurrently with any other
phase after P0. Its four classes currently **fail to encode the input at all** (0% pass), so every
fix is a direct pass-count win.

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-aromatic -b feature/roundtrip-aromatic release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-aromatic && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- prior diagnosis — `spec/handoffs/v0.3.6/S3-aromatic-perception.md`, `docs/KNOWN_LIMITATIONS.md`
  (porphyrinoid macrocycles "mostly resolved" — read what already round-trips before re-fixing)

### 3 · Verified code paths (all ENCODER-side)
- **kekulize_encode_crash (8):** `src/oinsmiles/utils/aromaticity.py:105` `kekulize_safe_sanitize`
  raises the two messages at `:141-143` ("cannot kekulize ... no quinoid ring to relax") and
  `:151-154` ("still invalid after de-aromatizing the quinoid ring(s) at atoms ..."). Helpers:
  `stuck_ring_atoms :30`, `clear_ring_aromaticity :59`, `dearomatize_stuck_rings :75`. Called from
  `xyz2mol.py:193,:234,:697`.
- **encode_crash_other (9):** `src/oinsmiles/utils/xyz2mol.py:625-634` raises "get_lig_mol failed for
  ligand fragment #{i} ...; cannot build TMC mol" when `get_lig_mol :426` → `_select_lig_mol :440`
  charge/carbene ladder returns `None` (`:625-626`). The `RuntimeError: Pre-condition Violation` is
  an uncaught native RDKit error on a perceived-but-invalid fragment.
- **macrocycle_perception (9, porphyrinoids):** rebuild loop `xyz2mol.py:1104-1128`; `AddBond` copies
  Kekulé types over an aromatic-flagged ring (hazard comment `:1116-1122`); re-perception `:931` can
  land on a different Kekulé than the input.
- **garbled_aromatic (7, `c=`):** `xyz2mol.py:893` `_repair_mixed_aromaticity` (called `:1165`)
  **early-returns unrepaired** when any bond is already AROMATIC-typed (`:909`), so a fragment mixing
  aromatic-typed and Kekulé-double ring bonds serializes as `c=`.

### 4 · Mission & scope guard
- **Re-triage first.** Some of these already round-trip on current `main` (porphyrinoid macrocycles
  are "mostly resolved" per `KNOWN_LIMITATIONS.md`; R3's `canonical_roundtrip_key` collapses
  aromatic/Kekulé). Confirm which goldens still crash on `c7edeeb6` before writing code — do not
  re-fix a solved case.
- **garbled_aromatic** is the cleanest target: narrow the `:909` early-return so a mixed
  aromatic/Kekulé fragment is actually repaired instead of bailed on.
- **kekulize/encode crashes**: extend the quinoid-relaxation / charge-carbene ladder to the failing
  skeletons. Watch: `SanitizeMol` rejects 4-coordinate neutral boron; `FastFindRings` ≠ SSSR.
- **macrocycle_perception**: any residual is a Kekulé-consistency issue in the rebuild loop
  (`:1104-1128`) — align the re-perceived Kekulé with the input for the specific porphyrinoids.
- Some `encode_crash_other` rows are **carborane** in disguise (fragment SMILES `'B'`) — those are
  wontfix; route to `docs`, don't chase.

### 5 · Owned files
`src/oinsmiles/utils/xyz2mol.py`, `src/oinsmiles/utils/aromaticity.py`. Fully disjoint from all other
phases.

### 6 · Gate
- Each still-crashing golden now encodes and round-trips (canonical key) via `--only` into
  `/tmp/rt-aromatic`. Use the **contract-mol simulation** (`Chem.Kekulize(tmc_mol)`, flags retained)
  for a ~1 s check against committed fixtures, per README.
- New guard tests under `tests/unit/`, each failing pre-fix.
- Full unit suite green (own baseline first); `ruff check` clean.
- Spot-check: a sample of currently-encoding aromatic complexes (ferrocene, fac-Ir(ppy)3, a
  metalloporphyrin that already passes) still encode byte-identically.

### 7 · Landing
Squash-merge into `release/v0.4.2` (see `SESSION_PROMPTS.md`).
