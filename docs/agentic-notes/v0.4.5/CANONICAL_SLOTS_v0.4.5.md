# Canonical coordination slots — v0.4.5 Lane 2

> **What this closes.** A donor's `{n}` slot integer used to be *which idealized template
> vertex its physical 3D direction vector best aligned to*, lex-**max**imized over a
> brute-forced rotation grid, with every remaining tie settled on the **raw input atom
> index**. That is a distortable, conformer- and file-order-dependent quantity, not a graph
> invariant. Behind `OIN_CANONICAL_SLOTS` (default **OFF**) the slot integers become a graph
> invariant, and the fragment order in the emitted string is re-derived from them.

Lever: `OIN_CANONICAL_SLOTS`. Default OFF. Levers-OFF output is byte-identical to the
v0.4.5 base, with **one accounted exception** — see §5.

---

## 1. The rule, and why it is correct

Fix the coordination geometry from the 3D fit. That determines the **isomer**: which
vertices are occupied and how they sit relative to each other (adjacent vs opposite =
cis/trans; facial vs meridional = fac/mer). Then choose, **among only the vertex labelings
the geometry's own symmetry permits**, the lexicographically minimal one.

The nuance that makes this correct rather than merely intuitive: **priority cannot be
applied independently of geometry.** Sorting ligands by priority and stamping them onto
slots 0,1,2,… would erase cis/trans and fac/mer, which are *defined* by which vertices
equal-priority ligands occupy. The relabeling must run inside the freedom the geometry
allows — the polyhedron's **proper-rotation group** (order 24 for an octahedron).

Proper rotations only. An improper operation maps a structure to its mirror image, so
folding over reflections would collapse enantiomers and flip the eta winding sense.
`derive_rotation_group` filters spanning vertex sets on `det > 0` for exactly that reason.
This is the guard rail the Y2 axial wave lacked when a tie-break sorted on a stereochemical
*sign* and silently made the token reflection-invariant.

## 2. Where it is applied, and why there

`utils/xyz2mol.get_oin_string`, as a post-pass on the **finished inline string**
(`canonical_slots.canonicalize_oin_slots`), not inside `_permute_and_serialize`'s lex-max
loop.

The inline string is exactly the representation `compare._parse_vertex_colors` already
reads. Reusing that function **verbatim** means the encoder's canonicalization and the
comparison key's canonicalization consume the same bytes through the same code, so they
agree *by construction* rather than by two parallel implementations staying in sync. It
also leaves the geometric fit, the eta RC1 content swap and the heading-atom tiers — all of
which read `item["slot"]` — running untouched beforehand.

Winding characters (`>` `<` `^`) are carried through **verbatim**. `_parse_vertex_colors`
folds `^` to `>` for *colouring* only; the emitted string keeps whatever the aligner
computed. Two independent reasons this cannot be invalidated by the relabel: winding is
measured against each ring's own metal→centroid axis (`oin_aligner._determine_winding`),
never against the slot's template direction; and the group contains only proper rotations,
which preserve circulation sense.

## 3. Fragment order

Re-derived from the canonical slots: `(minimum canonical slot, fragment text)`, with the
metal fragment pinned first (`fragments[0]` is a load-bearing project invariant — the
generator, the inline parser and the comparison key all assume it).

`minimum canonical slot` is by itself a **total** order on the coordinated fragments,
because a slot belongs to exactly one fragment; the text term only ever separates two
*uncoordinated* fragments. Slot-primary also matches the pre-existing convention (step 7 of
`get_oin_string` already sorted on `_sort_slot` primary), which keeps byte churn at
promotion time to the minimum.

This replaces `get_input_order_key` — the raw-XYZ-index tie-break — with a property of the
molecule instead of a property of the file. See §6 for what "replaces" means precisely.

*Deviation from the lane charter, on purpose.* The charter proposed
`(canonical body, min canonical slot, sorted canonical slots)`. Body-primary is equally
invariant, but it is not *more* invariant (min-slot is already a total order on the
coordinated fragments), and it reorders fragments relative to the pre-existing convention for
no gain — which would multiply the byte churn at promotion time for every molecule rather than
only the ones whose slots actually move. `sorted canonical slots` is likewise unreachable: two
fragments cannot share a minimum slot.

## 4. The three order-dependent seams that had to be neutralized together

Slot relabeling alone is not enough, because two upstream steps make the *frame* itself a
function of the input numbering, and a different frame produces a different geometric fit:

| seam | was | now (lever ON) |
|---|---|---|
| `_align_to_pai` pivot | `np.min(candidates)` — lowest **file** index among atoms tied for max distance from the metal | `_canonical_pivot`: `(-mass, sorted multiset of interatomic distances)` — an atomic property plus a rigid-motion invariant |
| `_align_to_pai` Z-sign | `sum(z_i * (i+1)**3)` — a weight built from the atom's **position in the file** | `_canonical_z_sign`: mass-weighted **odd** moments in z, tried in order; plain sums, so unchanged by renumbering, and odd, so they actually decide between `+z` and `-z` |
| `get_input_order_key` | `min(original atom indices)` | donor position inside the fragment's own canonical SMILES (`s_idx`, off `_smilesAtomOutputOrder`) |

## 5. The rotation group is now unified — and this is the one unconditional change

`oin_aligner._brute_force_symmetries` used to build the group by brute-forcing Euler
triples from the fixed grid `[0, 90, 120, 180, 240, 270]`. **That grid cannot express a
72° five-fold rotation**, so on PBP it found **2 of the 10** proper rotations and the
encoder could not canonicalize a pentagonal-bipyramidal equatorial labeling at all. It
agreed with `derive_rotation_group` on the other ten geometries and never invented a
non-rotation.

It now delegates to `canonical_slots.derive_rotation_group` (memoized per template), so the
encoder, the comparison key and the canonical-slot post-pass read **one** vertex table and
**one** group derivation. That closes open debt TD-005.

`tests/unit/test_canonical_slots.py::TestAlignerSymmetriesAreUnified` pins both halves:
agreement on every geometry, and the specific C5 generator `(0,1,3,4,5,6,2)` present in the
PBP group so an angle-grid search cannot silently come back.

⚠ **This delegation is NOT behind the lever**, so on PBP complexes it changes the
default-path bytes: eight extra permutations enter `_permute_and_serialize`'s lex-max.
Measured incidence is in §7. Everywhere else it is provably a no-op.

## 6. What "delete the `get_input_order_key` tie-break" actually means

Deleting it outright would change levers-OFF bytes, which the release's default-OFF
discipline forbids. It is therefore *superseded* rather than textually removed, and the
supersession is total on the canonical path:

- Its **second** use (step 7, the final fragment sort) is overwritten wholesale — the
  post-pass re-derives fragment order from the canonical slot integers.
- Its **first** use (step 5, before the aligner runs) only ever separates fragments that
  already tied on mass, binder mass **and** body SMILES — i.e. two copies of one ligand.
  Two fragments with identical bodies and identical donor elements get identical vertex
  *colours*, so swapping them leaves `vcolor` unchanged and the lex-min labeling identical.
- Under the lever its value is no longer a file index at all.

The empirical statement of the same claim is the acceptance test: permute the input atom
order and the emitted string must be byte-identical (`tests/unit/test_canonical_slot_invariance.py`).

## 7. Measured

Instrument: `tools/canonicality_probe.py --n 150 --trials 2` (seed 42, so every arm samples
the same molecules), generator-free, on branch `swimlane/v045-lane2`. Transforms hold the
molecular graph fixed, so the correct answer is byte-**identical**, not merely similar.

Lane 2 declares a hard dependency on Lane 1 (the canonical ligand body *is* the colour the
slot lex-min minimizes over), so the arm that matters is the stacked one. Lane 1's levers
are `OIN_CANONICAL_BODY` + `OIN_CANONICAL_PERCEPTION`; `OIN_STABLE_METAL_AC` is this
branch's own AC fix.

| arm | byte-stable | drifted | `rdkit_canonical` | `slot_renumber` | `fragment_reorder` | `winding_star_drift` | `rotate` drift | key-level defects |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all levers OFF | 88/149 (59.06%) | 61 | 39 | 25 | **0** | **0** | **0** | 23 |
| `OIN_CANONICAL_SLOTS` only | 94/149 (63.09%) | 55 | 28 | 29 | **0** | **0** | **0** | 23 |
| Lane 1 levers only | 95/150 (63.33%) | 55 | 34 | 23 | **0** | **0** | **0** | 18 |
| **stacked (Lane 1 + slots)** | **100/150 (66.67%)** | 50 | **19** | 32 | **0** | **0** | **0** | **18** |

(The two OFF-side arms have 149 encoded rather than 150 because `OIN_STABLE_METAL_AC` also
recovers one input the base encoder fails on. Read the byte-stable *counts*, 88 → 100.)

Per-molecule attribution, not just counters:

| comparison | fixed | **broke** |
|---|---:|---:|
| OFF → slots only | 6 | **0** |
| Lane 1 only → stacked | 5 | **0** |

**`slot_renumber` rising is not a regression — it is reclassification.** The probe's
taxonomy tests `slot_renumber` *before* `rdkit_canonical`, so a molecule whose body drift
closes while its slot drift remains *moves into* the finer class. Measured: 8 molecules
`rdkit_canonical → slot_renumber` on the OFF→slots step and 11 more on the
Lane1→stacked step, with zero molecules going the other way and zero going stable→drift.

**The slot lever moves nothing at isomer level, in either direction** (key-level defects
23 → 23 alone, 18 → 18 stacked). That is the over-folding guard at corpus scale: had the
canonicalization merged a real isomer, the comparison key would have registered it.

### Exactly which pinned goldens the lever moves — and proof each is a relabeling

With `OIN_CANONICAL_SLOTS` alone, three of the six pinned fixtures move and three are already
canonical. **Every moved string has an unchanged `canonical_roundtrip_key`**, i.e. the lever
relabeled the isomer without becoming a different one.

| fixture | off | on |
|---|---|---|
| CisPlatin | `[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}` | `[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}` |
| TransPlatin | `[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}` | `[Pt_SPL].N{0}.[Cl]{1}.N{2}.[Cl]{3}` |
| fac-Ir(ppy)₃ | `…c{0}…n{3}1.c{5}…n{1}1.c{2}…n{4}1` | `…c{0}…n{5}1.c{2}…n{1}1.c{4}…n{3}1` |
| mer-Ir(ppy)₃, Ferrocene, Cis-PtCl₂(en) | — | unchanged |

Cisplatin is the clearest case. The old labeling put the chlorides on slots 0,1 only because
the fragment sort ranks by descending mass. The lex-min runs over vertex *colours* and
`"N" < "[Cl]"` bytewise, so the canonical labeling puts the amines on the lowest vertices.
**Cis is still cis** — the chlorides land on 2 and 3, which are adjacent.

All of this is pinned by `TestLeverOnGoldens` (moved strings *and* the key-invariance
assertion *and* the unmoved set, so gratuitous churn is caught too), which turns the
promotion A/B into a diff against a committed expectation rather than a judgement call.

Consequence to expect at promotion time: `tests/unit/test_regression_stability.py` pins the
**default path** and does not clear the environment, so running the suite with
`OIN_CANONICAL_SLOTS` exported fails exactly those 3 goldens (`test_cisplatin`,
`test_transplatin`, `test_fac_irppy3`) and nothing else. Left as-is on purpose: that file is
the project's default-path pin and hardening it is a change to a shared guard, not this
lane's call. `TestLeverOffIsByteIdentical` re-checks two of the same goldens with the
environment force-cleared, so golden coverage survives either decision.

### Rotation-group unification: measured blast radius

Deterministic encoder-vs-encoder byte diff over the same 150 molecules, derived group vs the
historical Euler grid, levers OFF:

- geometry incidence: `SPL` 46, `OCT` 33, `TET` 32, `LIN` 11, `TPL` 9, `TBP` 7, `SPY` 5,
  `TPY` 4, **`PBP` 2**;
- **byte-different: 1/149 (0.67%)** — `IMACOO_comp_0 [Zr_PBP]` only. Every other geometry is
  a no-op, as the docstring claims.

**It is left unconditional, deliberately.** `compare.py`'s `_polyhedron_signature` has used
the *derived* group since v0.4.4 and was never gated, so on the default path the comparison
key already enumerated all 10 PBP rotations while the encoder enumerated 2. That is a live
encoder/key disagreement, and unifying removes it rather than introducing anything. The
change folds over proper rotations only, so no isomer information can be lost, and every
pinned golden stays byte-identical. If the release wants strict corpus-wide byte-identity
with levers off, this is the one line to gate — flagging it rather than deciding it silently.

## 7a. WHERE THE PLAN WAS WRONG — the residual class, and why the post-pass cannot reach it

> **The plan's premise was that "the permutation `_polyhedron_signature` achieved its lex-min
> with **is** the canonical slot relabeling." That is necessary but NOT sufficient, and the
> acceptance target `slot_renumber → ~0` is therefore not reachable at this seam.** Measured,
> not argued.

`slot_renumber` residual in the stacked arm is **32 molecules**, and the mechanism is
uniform: **32/32 are `same_vcolor_identical`** — `compare._parse_vertex_colors` returns the
*identical* `{slot: colour}` map for the base string and the drifted one.

That is fatal for a post-pass that derives its relabeling from `vcolor` alone: identical
input, identical permutation, so applying it to two differently-labeled strings **preserves
the difference exactly**. Not a bug in the relabeling — a limit of its input.

The reason the colours are identical is a deliberate design choice one layer down.
`_parse_vertex_colors` colours *every* donor of a ligand with that ligand's **whole** body,
with no chelate grouping, precisely so that "a swap between two same-coloured donors is
invisible" (its own docstring). That blindness is what lets true conformers collapse in the
comparison key while fac/mer stay distinct. It is correct for a key. It is **insufficient for
an encoder**, which has to emit one string and so must decide something the key is entitled
not to care about: *which donor atom of a chelate holds which slot integer.*

Worked example, `AGUKOD_comp_0` (Rh, SPL, a COD ligand spanning two cis vertices):

```
base: [Rh_SPL].[CH]{1}1=[CH]{1>}CC[CH]{0>}=[CH]{0}CC1.…{2}….[Cl]{3}
got : [Rh_SPL].[CH]{0}1=[CH]{0>}CC[CH]{1>}=[CH]{1}CC1.…{2}….[Cl]{3}
```

Bodies byte-identical; `vcolor` identical (`{0: (cod, C, ">"), 1: (cod, C, ">")}`); only which
alkene arm carries which integer moved. And the transposition `(0 1)` that relates them is
**not an element of SPL's rotation group** (D4 contains `(1 3)`, `(0 2)`, `(0 1)(2 3)`, … but
not `(0 1)` alone), so the group-theoretic fold cannot see it either.

### Atom-level verdict: 23 benign, 7 that are worse than a canonicality defect

`tools/slot_drift_mechanism.py` resolves each residual pair to the donor atoms involved and
asks whether they are interchangeable, using `CanonicalRankAtoms(breakTies=False)` symmetry
classes (`breakTies=True` settles ties between symmetry-equivalent atoms on the *input index*
— the very dependence under measurement — so it is the wrong instrument; `includeChirality`
defaults True, so two constitutionally-equivalent branches with different configurations
correctly land in different classes):

| verdict | n | meaning |
|---|---:|---|
| `automorphism` | **23** | the two donor atoms are in one symmetry class of their own fragment, so both labelings denote the **same molecule**. The encoder simply lacks a deterministic choice. |
| `DISTINCT_donors` | **7** | the slots land on **inequivalent** donor atoms. One of the two strings is therefore **wrong** — a soundness defect, not a canonicality one. |
| unparsable | 2 | fragment would not sanitize (borane-cluster class). |

The 7: `IMACOO_comp_0`, `RUBTIS_comp_0`, `VEXHIR_comp_0`, `ZACFER_comp_0`, `HAVGIW_comp_0`,
`KISQAG_comp_0`, `ZOSNUS_comp_0`.

### Root cause of the 23, named

The physical donor→vertex map comes from the 3D fit and is order-free. What moves is **which
of two symmetry-equivalent donor atoms RDKit writes first in the canonical SMILES** —
`_smilesAtomOutputOrder`, i.e. `CanonicalRankAtoms(breakTies=True)`, whose tie between
symmetry-equivalent atoms is settled on the input index. The canonical *string* is invariant;
the *map from string position to atom* is not. So the same physical assignment gets written
with the two integers exchanged.

### The fix, and why it was NOT taken here

Fold, in addition to the rotation group, the permutations that exchange donors **within one
fragment** that are (a) in the same `breakTies=False` symmetry class and (b) the same colour.
Scoped that narrowly it does not touch metal Δ/Λ (which lives in the arrangement *across*
fragments, still folded only by proper rotations) and it does not touch the 7
`DISTINCT_donors` cases (different classes → never folded).

It is not implemented because **over-folding is this lane's entire risk** and this widens the
fold beyond the geometry's own symmetry — the one boundary the whole design rests on. It
needs its own guards (a Δ/Λ tris-chelate divergence fixture at minimum) and a product call,
not a unilateral extension by the agent that found it. Recorded here so it can be picked up
in an hour rather than re-derived.

**Consequence for the acceptance criterion.** `slot_renumber → ~0` is not achievable by
canonicalizing slot *labels*; it needs the *atom→label* assignment made order-invariant too.
What this lane does deliver, measured: **+12 byte-stable molecules (88 → 100 of ~150), 0
regressions, 0 isomer-level change, `fragment_reorder` and `winding_star_drift` held at 0.**

## 8. The export for Lanes 5 and 6

```python
from oinsmiles.oin.canonical_slots import canonical_slot_map

canonical_slot = canonical_slot_map(oin_string)[slot_as_it_appears_in_that_string]
```

`oin_string` is what `get_oin_string` / `XYZToSMILES.convert` returned. The map is the
**identity** when `OIN_CANONICAL_SLOTS` is on, because the encoder has already applied it —
which is the point: a caller written against this helper is correct with the lever either
way and stays correct after promotion. It is idempotent on an already-canonical string.

**Do not use `canonical_slot_permutation(geo, vcolor)` for this question.** The lex-min
*signature* is invariant, but the permutation achieving it need not be unique: when the
coloured polyhedron has a nontrivial colour-preserving rotation stabilizer several
permutations tie, and that function breaks the tie on the permutation tuple, which is a
property of the *incoming* labeling. `canonical_slot_map` refines the tie on the rendered
output, which is invariant by construction (proof in its docstring: for `L' = g·L` the
candidate sets are equal because `p ↦ p·g` is a bijection of the group).

**What is and is not unique.** The canonical *string* is unique, always. The per-donor map
is unique only up to the coloured polyhedron's **rotational** automorphism group — and
exactly up to it. fac-M(ppy)₃ has a real C3 axis through its three equivalent ligands, so
three donor labelings are equally canonical while the emitted string stays single-valued.
That is not a defect to route around:

- every member of that group is a **proper** rotation, so a descriptor that flips under
  reflection — metal Δ/Λ, an eta winding, an axial sign — takes the *same* value on all of
  them;
- a descriptor that is *not* invariant under the automorphism is not a property of the
  molecule at all; it is a property of which interchangeable ligand you chose to call
  "first".

Derive from the canonical arrangement and let this map join a donor to it.

## 9. Known limits / not addressed here

- **The metal's own `@SP1`/`@OH1` tag would be a parity against a donor ordering the
  post-pass renumbers** — but **measured 0 of 150 molecules emit a metal `@` tag at all**, so
  there is nothing stale to break today. (`normalize_oin_for_comparison` strips metal `@`
  stereo anyway, so no current measurement could see it either.) Making the metal descriptor
  exist and be reproducible is Lane 5's charter, and Lane 5 must derive it from the canonical
  arrangement (§8) rather than from pre-relabel slots.
- Unknown geometry tags fold with the identity only. Conservative by design: an unknown
  geometry can *over*-split (miss a benign rotation) but can never wrongly merge two
  isomers, so it cannot reintroduce fac/mer blindness.
- Bond-order / aromaticity perception order-dependence (`CEBVIR_comp_0`) and the chiral-tag
  parity class (`FEQFIS_comp_0`) are **not** slot problems and survive this lever; they are
  Lane 1's and Lane 8's respectively. Confirmed in
  `docs/agentic-notes/v0.4.5/RENUMBERING_INSTABILITY_v0.4.5.md`.
- `DUDREA_comp_0` is closed by `OIN_STABLE_METAL_AC` on this same branch, and it was never
  a slot problem either — the AC valence-capping loop iterated in atom-index order. Do not
  reopen it.
