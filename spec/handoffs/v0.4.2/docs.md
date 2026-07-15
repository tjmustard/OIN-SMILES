# ▶ START HERE — docs phase (v0.4.2 round-trip accuracy wave)

**Launch a fresh Claude Code session in the main checkout and hand it this file.** A lightweight,
file-disjoint phase — owns `docs/KNOWN_LIMITATIONS.md` + `spec/handoffs/v0.4.2/wontfix-carboranes.md`.
Runs concurrently with anything after P0; lands whenever ready. It **collects** the notation/harness
residuals the fix phases route here.

### 1 · Create and enter your worktree
```bash
git -C /home/tjmustard/Documents/GitHub/OIN-SMILES worktree add \
  /home/tjmustard/Documents/GitHub/OIN-SMILES-docs -b feature/roundtrip-docs release/v0.4.2
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-docs && uv sync
```

### 2 · Read these (main checkout)
- shared protocol — `spec/handoffs/v0.4.2/README.md`; floor — `spec/handoffs/v0.4.2/BASELINE.md`
- current — `docs/KNOWN_LIMITATIONS.md` (already documents porphyrinoids, carboranes, atom-stereo
  residuals, no_conformers), `spec/handoffs/v0.3.6/wontfix-carboranes.md`

### 3 · Mission — document, don't fix
Add/refresh `docs/KNOWN_LIMITATIONS.md` entries for the classes v0.4.2 deliberately does not fix, and
the residuals fix phases hand off:
- **carborane_unsupported (79):** ENCODER wontfix — `get_lig_mol` (`xyz2mol.py:625-634`) → `AC2mol`
  is a two-centre bond-order solver; 3c2e boron cages need multi-centre bonding it cannot represent
  (`KNOWN_LIMITATIONS.md:144-150` already states this; confirm the count and add the current member
  list). Mirror `spec/handoffs/v0.3.6/wontfix-carboranes.md` into `spec/handoffs/v0.4.2/wontfix-carboranes.md`.
- **donor_H nitride/ammine ambiguity (S1 residual):** the +3 `N{n}` rows where nitride `[N]` and
  ammine `[NH3]` both serialize to `N{n}` and "ammine reading wins" (`metallogen_adapter.py:214-219`)
  — a notation limit needing an OIN format change. Document the residual set S1 hands over.
- **high_rmsd FF floor + `--quick` timeout (S7 residual):** the 1.0 Å gate runs on chemically-correct
  FF-only geometries; `--quick` timeouts are a 30 s/max_attempts=10 budget artifact. Document that
  these are *not* accuracy failures and what a full-fidelity run recovers.
- **UncoordinatedFragmentError (S7 residual):** outer-sphere counterions/solvents are not
  representable in MetalloGen m-SMILES — a representation limit.
- **Irreducible generator stereochemistry:** QOFTOU-class (builds rac/meso non-deterministically) and
  any winding row S5 could not make deterministic.

### 4 · Gate & landing
- Each residual class has a `KNOWN_LIMITATIONS.md` entry naming the layer that would have to change,
  with a current member list (pull from a **fresh** `classify_failures.py` run on a **copy**).
- No source changed. Squash-merge into `release/v0.4.2` (see `SESSION_PROMPTS.md`).
