# v0.4.6 — H-faithful canonical body, and two refuted hypotheses

Branch `swimlane/v046-hfaithful`, off tagged v0.4.5. **Nothing here is merged.**

The v0.4.5 close-out ranked four items toward the accuracy goal. Two of the four turned out to be
wrong as stated, and finding that out is the main result of this branch.

## 1 · H-faithful canonical body — CORRECT, but NO measured benefit

**The mechanism was real.** `canonical_body_emit` had two writes, both plain `Chem.MolToSmiles`:
the intermediate that feeds the reparse, and — the one that matters — the **final emit whose
output becomes the body**. So `xyz2mol.py:1710` computed an H-faithful string and `:1736`
overwrote it. That is exactly how `OIN_CANONICAL_BODY` "undid" `OIN_H_FAITHFUL`, as recorded in
`levers.py::_HELD_OFF`. Both writes now route through `h_faithful_smiles`.

**The benefit was not.** A/B over the 45-molecule `Atom count mismatch` population, comparing the
OIN-implied atom count against the input XYZ (generator-free):

| arm | match | mismatch |
|---|---|---|
| `OIN_H_FAITHFUL=0`, canonical body ON | 8 | 37 |
| `OIN_H_FAITHFUL=1`, canonical body ON | **8** | **37** |

Identical. **The atom-count class is not a serialization defect.** `h_faithful_smiles` guarantees
only that a string re-reads with the hydrogens it was *written* with; if the perceived molecule
already disagrees with the input XYZ, the string faithfully describes a wrong graph and no
emit-side change can help. This is the phantom-hydrogen class, and it lives in perception
(`get_lig_mol` / `AC2BO`), not in `MolToSmiles`.

**Why the change is kept anyway:** it removes a genuine inconsistency — one lever silently
undoing another — and it is a prerequisite for `OIN_H_FAITHFUL` ever mattering. It is **not** an
accuracy win and must not be counted as one.

**Verified:** byte-identical on all 61 fixtures against unmodified `main` (`#DONE 61` sentinel on
both arms, outputs diff-clean), and the full suite is 838 tests OK / 3 skipped / 4 xfail.

## 2 · Restoring the locked-donor tag through the reparse — REFUTED

Four lines. The correspondence is already available: `_reparse_once`'s Guard 2 proves
`donors[k] <-> new_donors[k]` (same element, same heavy degree). It makes P3 emit under
`OIN_CANONICAL_BODY`, and **POJJOP passes**.

It is still wrong. Setting a chiral tag *after* the sanitize introduces a stereocentre the
canonical ranker did not account for, which moves the canonical **write order** — and `@`/`@@` is
a parity relative to that order, not an absolute label. On `RIFGUJ_comp_2` (three Cu-bound amines
on one cyclohexane) the three ring-**carbon** tags then flip between a structure and its mirror.

The geometry says they must not. `AssignStereochemistryFrom3D` + `rdCIPLabeler` label those
carbons lowercase **`s`** — pseudo-asymmetric, a *relative* (all-cis) descriptor — and they read
`s` identically for the structure **and** its reflection:

```
atom   6  base=s  mirror=s   same
atom   8  base=s  mirror=s   same
atom  10  base=s  mirror=s   same
```

So the restoration silently rewrote stereochemistry that must not move. The lane's multi-centre
mirror guards caught it; single-centre POJJOP could not — the Y2 lesson intact. **v0.4.5's
decision to defer this rather than rush it was therefore correct.**

Reverted, with the mechanism written into `canonical_body.py::_reparse_once`,
`levers.py::_HELD_OFF` and the test module, plus a new lever-independent guard
`test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric` so the four-line "obvious fix"
cannot be silently re-attempted.

A correct fix must preserve the tag **without perturbing the ranking**: keep the donor bracketed
through the sanitize, or re-derive parity from the parent geometry once the write order is fixed.

## 3 · The timeout bucket masks defects — REFUTES "the gap is mostly compute"

24 molecules (seed 42, stratified 12 `UFF_1` + 12 `g-xTB_1`) that all hit the 300 s wall in full
mode, re-run on the cheaper `--quick` path:

| outcome | n |
|---|---|
| SUCCESS | 6 (25%) |
| String mismatch | 6 (25%) |
| Atom count mismatch | 6 (25%) |
| MetalloGen failed | 6 (25%) |
| **timed out again** | **0** |

Only ~25% is compute-limited. More compute buys ~44 of 936 molecules (~4.7%), not the 174 the
timeout count implies.

⚠ The probe's `elapsed_s` values are unusable — it ran alongside the 6-shard 5k sweep at load
21–33, and `tools/v045_state.sh` says wall-clock is meaningless above ~12. Only the pass/fail
outcomes survive, and only because none of the 18 failures was a timeout.

## Revised ranking (the previous one is superseded)

1. **`OIN_BORON_CAGE`** — the only item with a measured accuracy gain: 0/36 → **34/36** encodes on
   the boron population, 0.2–4.2 s each. Costs 14 silent false passes becoming honest failures.
   Needs a full sweep, not a mid-release flip.
2. **Perception-side hydrogen** — the atom-count class, now correctly located in `get_lig_mol` /
   `AC2BO` rather than in serialization. Unscoped.
3. **MetalloGen generation failures** (~80 molecules). Unscoped.
4. **String mismatch** (~55 including timeout-masked). Unscoped.
5. **Compute** for the ~44 genuinely timeout-limited. Buys the least of any option.

## On the 100% target

Not reachable from here, and not only for effort reasons. `xyz2AC_obabel` can perceive a
genuinely **different graph** between two conformers near the covalent-radius + 0.45 Å cutoff
(`xyz2mol_local.py`), which bounds what any canonicalization can achieve — this is stated as
out-of-scope in the v0.4.5 plan and remains true. The <30 s target is much closer: already 83.6%
of succeeding molecules, median 6.9 s, p90 50 s.
