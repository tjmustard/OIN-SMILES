# v0.4.2 BASELINE — clean single-commit floor (PLACEHOLDER)

> **Not yet produced.** This file is written by **P0** (`P0-baseline.md`). Until P0 lands it, the
> only failure data is the mixed-provenance, `--quick` accumulator in
> `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0/` — **not a valid floor** (do not quote its headline %).

When P0 lands, this file will contain, all on the single baseline commit **`c7edeeb6` (= tag
`v0.4.1`)**:

- The exact sweep command, rdkit version, date.
- **Per fixable class** (donor_H fixable subset, H_on_terminal_oxo_imido, geometry_NON,
  geometry_or_fragment_change, winding_flip, EZ_bond_stereo, atom_stereo, encode_crash_other,
  kekulize_encode_crash, macrocycle_perception, garbled_aromatic, `[S@SP3]` subset): current count +
  the frozen **goldens** (4–8 molecules each, reproduced on `c7edeeb6`) that each phase A/Bs against.
- The **must-not-regress passing set** — molecule IDs stamped `commit_id == c7edeeb6` **and**
  `status == success` — and, separately, the **untrusted-provenance** (`5538b722-dirty`) passers.
- The artifact-class counts (timeout / high_rmsd / carborane / no_conformers) as **context only**,
  marked "not part of the floor; S7/docs own them."
- **No headline percentage.** The floor is a *set of molecule IDs*, not a number.
