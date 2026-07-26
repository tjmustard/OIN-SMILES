# Encoder instability under input atom renumbering — measured 2026-07-25

> **Headline: 13.0% of sampled molecules emit a different absolute stereochemistry when the
> same 3D structure is presented with its atoms in a different order.** This is a soundness
> defect, not a formatting one, and it was not on the v0.4.5 plan.

Measured with `tools/canonicality_probe.py` at `main` @ `20044883`, all levers OFF.
Sample: 225 molecules (seed 42) from the 25,197-basename tmCAT-tmPHOTO corpus, 2 trials per
transform. 223 encoded successfully.

## What was done

The probe takes one input structure and re-presents it three ways, **holding the molecular
graph fixed** so the expected answer is byte-identical, not merely similar:

| transform | what changes |
|---|---|
| `rotate` | a random **proper** rotation (det = +1) of the coordinates |
| `renumber` | the order atoms appear in the XYZ file |
| `both` | renumber, then rotate |

Rotations are forced proper because an improper operation mirrors the structure, which
legitimately changes a chiral molecule's encoding.

## Results

| | count | share |
|---|---:|---:|
| byte-stable across all transforms | 125 / 223 | **56.1%** |
| drifted | 98 / 223 | **43.9%** |
| **of which: comparison KEY also changed** | 47 / 223 | **21.1%** |

Drift by transform: `renumber` 132, `both` 123, **`rotate` 0**.

### `rotate` produces exactly zero drift

Across 225 molecules × 2 trials, not one structure changed its OIN string under a random
proper rotation. **The encoder is fully orientation-invariant.** `_align_to_pai` does its job
for orientation. Every defect below is an *atom-numbering* dependence.

### Severity classes (under pure renumbering)

| class | count | share |
|---|---:|---:|
| **stereo-only flip** — string identical after deleting `@`/`@@`, tags differ | **29** | **13.0%** |
| other drift — skeleton, slot numbering, or aromaticity also differs | 79 | 35.4% |
| geometry classification changed (`[M_XXX]`) | 1 | 0.4% |
| aromaticity perception collapsed (>2 atoms) | 8 | 3.6% |

The stereo-only class is measured **strictly**: the two strings are byte-identical once
chiral tags are removed. Because SMILES chirality is defined relative to the order neighbours
appear *in the string*, and that order is identical here, the differing tags denote genuinely
different configurations. At least one of the two is therefore wrong.

## Worked examples

**`FEQFIS_comp_0`** — verified the two inputs are a pure permutation (identical sorted
coordinate multiset, so no geometry change and no mirroring):

```
BASE : [Au_LIN].C[C@@H](c1ccccc1)N([C@@H](C)c1ccccc1)p{0}1oc2c(...)...[Cl]{1}
RENUM: [Au_LIN].C[C@H](c1ccccc1)N([C@@H](C)c1ccccc1)p{0}1oc2c(...)...[Cl]{1}
                    ^^^^^
```
Byte-identical apart from one tag.

**`DUDREA_comp_0`** — the **coordination geometry classification changed**, and a different
set of atoms is recorded as coordinated:
```
BASE : [Y_SPY]. ... .[BH3]{2}.[H]{3}.[BH3]{4}.[H]
RENUM: [Y_TET]. ... .[BH3]{2}.[BH3]{3}.[H].[H]
```
Most likely route: the `(i+1)**3` Z-moment weighting at `xyz2mol.py:971` flips Y/Z under
renumbering, and the geometric template fit then selects a different polyhedron. The v0.4.5
plan called this seam "latent"; it is not latent.

**`CEBVIR_comp_0`** — aromaticity perception collapses entirely: `c1c(F)c(F)...n{0}c(...)`
becomes `C1C(F)C(F)...[CH]...`. Consistent with `AC2BO`'s "arbitrary resonance form" and
`get_UA_pairs`' non-unique `nx.max_weight_matching` being atom-order dependent.

## Why the Y1/Y2/Y3 audit did not find this

The injectivity audit asked **"does the encoder *separate* two enantiomers?"** — mirror-twin
collision probes. It never asked **"does the encoder *consistently report* one enantiomer?"**
Those are different questions, and the second was never instrumented. A mirror-twin probe
compares two *different* structures; this probe compares one structure with itself.

## Why the round-trip sweep did not find this either

The sweep compares an input encoding against a *generated* structure's encoding. The
generator builds coordinates from the OIN string, so both sides inherit whatever
configuration the encoder chose on that particular input ordering — the error is common-mode
and cancels. Only re-presenting the *same* input differently exposes it.

## Consequence for historical numbers

The comparison key is the harness's acceptance predicate. It changed under renumbering for
**21.1%** of sampled molecules, so reported round-trip accuracy carries a systematic error
term that has never been accounted for. This does not invalidate the relative A/B results
(both arms share the input ordering), but absolute accuracy figures should be read with it in
mind.

## Scope call

This is squarely within "the OIN string is a faithful 1D hash" and is more serious than the
canonicality drift v0.4.5 set out to fix, so it is being investigated as an added lane rather
than deferred. Suspects, in order:

1. `AC2BO` / `get_UA_pairs` order-dependence changing perceived bond orders, hence CIP
   priorities (`utils/xyz2mol_local.py:800`, `:542`) — shared root with Lane 1 step 2.
2. `core/chirality.py` `CIPAssigner` / `ChiralityRecoveryUtility` ordering assumptions.
3. `_align_to_pai`'s index-dependent pivot and `(i+1)**3` Z-sign (`utils/xyz2mol.py:941`,
   `:971`) — confirmed live by `DUDREA_comp_0`.

## Follow-up, measured 2026-07-25 (later the same day)

### One of the three mechanisms is now diagnosed to the line and fixed

**`DUDREA_comp_0` — CLOSED behind `OIN_STABLE_METAL_AC`** (`swimlane/v045-lane2` @ `8bf9df61`).
It was never a slot-labelling problem, which is why Lane 2's slot post-pass could not touch it.

`xyz2AC_obabel`'s distance pass (`utils/xyz2mol_local.py:1194-1202`) is order-**free**: a symmetric
comparison of the distance matrix against the covalent-radius sum. The **only** order-dependent
step in AC perception is the valence-capping loop that follows, which iterated
`for i in range(num_atoms)` — in input atom order. Capping atom *i* removes a bond, lowering some
atom *j*'s count, so whether *j* still needs capping depends on whether *i* was visited first.

This molecule is a Y borohydride whose bridging hydride is bonded to **both** B and Y, exceeding
H's valence of 1. Cap Y first → the Y–H bond survives → metal degree 5 → geometry `SPY`. Cap that
H first → it drops Y–H instead → degree 4 → `TET`. The losing contact was **not** the shortest:
the metal has six hydrides at 2.298 / 2.300 / 2.328 / 2.379 / 2.408 / 2.421 Å and the one that
flipped was 2.328 Å — so this was never "shortest wins", it was iteration order.

Measured: the perceived AC differed in **3 of 8** random renumberings with the lever off and
**0 of 8** with it on, metal degree settling consistently at 7.

⚠ Not promotable on this evidence. It changes **perception**, and capping the metal first can only
*keep* bonds the old order discarded — so coordination numbers can rise and geometries reclassify
elsewhere in the corpus. Needs a corpus A/B before any default flip.

### A hypothesis recorded earlier in this document's programme was WRONG

I had reasoned, and told the Lane 8 agent in writing, that the 13% stereo-flip class was probably
downstream of order-dependent bond-order perception and so would be fixed by Lane 1's
`OIN_CANONICAL_PERCEPTION`. **Measured on the three worked examples, 3 trials × 3 transforms:**

| levers | byte-stable | key-level defects remaining |
|---|---|---|
| all OFF | 0/3 | **3** |
| `OIN_STABLE_METAL_AC` | 0/3 | **2** |
| `+ OIN_CANONICAL_PERCEPTION` | 0/3 | 2 |
| `+ OIN_CANONICAL_BODY + OIN_CANONICAL_SLOTS` (all four) | 0/3 | **2** |

Only the metal-AC fix closes anything. **`FEQFIS_comp_0`'s stereo flip and `CEBVIR_comp_0`'s
aromaticity collapse survive all four levers.** So the stereo class does **not** share a root with
Lane 1's perception work, and Lane 8's independent investigation is genuinely required rather than
redundant. Its WIP (`swimlane/v045-lane8` @ `8fdccb55`) should be finished, not folded into Lane 1.

Scope caveat on this table: n = 3, deliberately the **hardest** hand-picked cases. It does not
refute Lane 1's general result (6 fixed / 0 regressed over 250 molecules); it refutes only the
specific claim that Lane 1's lever closes the stereo-flip class.

## Reproducing

```bash
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py --n 300 --trials 2 \
    --out tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-canonicality/baseline-main
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py \
    --only FEQFIS_comp_0,DUDREA_comp_0,CEBVIR_comp_0 --trials 3 -v
```
