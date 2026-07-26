# Wave B — canonical coordination slots (Lane 2), the release's one real dependency

**Purpose of the wave:** make a donor's `{n}` slot integer a **graph invariant** instead of a
property of how the input XYZ happened to be oriented and numbered — by choosing, among only the
vertex labelings the coordination geometry's own symmetry permits, the lexicographically minimal
one. Wave B is a wave of exactly one lane because Lane 2 **hard-depends** on Lane 1: the canonical
ligand body *is* the vertex colour that the slot lex-min minimizes over.

---

## ELI5

In an OIN string each atom that binds the metal carries a small integer — its "slot" — saying which
corner of the coordination polyhedron it sits on. Those integers are what make `cis` different from
`trans` and `fac` different from `mer`, so they cannot simply be sorted into a tidy order; that
would erase the very isomer information they exist to record. Before v0.4.5 the integer came from
fitting the real 3D directions onto an idealized template and then breaking every remaining tie on
"whichever atom appeared first in the file" — so two orientations of the same molecule got two
different but equally correct labelings, and the strings differed. The fix is to keep the geometry's
answer about *which corners are occupied* and only shuffle labels within the freedom the polyhedron's
own rotational symmetry allows, then pick the alphabetically smallest result. Two safety rules make
this sound: only **proper** rotations are allowed (a reflection would turn a molecule into its mirror
image and quietly merge two enantiomers), and the shuffle is applied as a final pass on the finished
string rather than surgically inside the code that built it.

## The wave, visually

```
   ┌──────────────────────────── WAVE A ────────────────────────────┐
   │  LANE 1  canonical ligand body  ──── the VERTEX COLOUR ────────┼──┐
   │          (OIN_CANONICAL_BODY, OIN_CANONICAL_PERCEPTION)        │  │
   └────────────────────────────────────────────────────────────────┘  │
                                                                 HARD │ DEPENDENCY
                                                                      ▼
  ╔══════════════════════════ WAVE B — LANE 2 ═══════════════════════════════════╗
  ║                                                                              ║
  ║  STEP 0 (shipped FIRST, deliberately RED on one geometry)                     ║
  ║    tests/unit/test_canonical_slots.py::TestAlignerSymmetriesAreUnified         ║
  ║      derive_rotation_group  vs  _brute_force_symmetries Euler grid            ║
  ║      ✔ agrees 10/11 geometries        ✖ PBP: 2 of 10 proper rotations         ║
  ║      ⇒ the PBP fix becomes a MEASURED improvement, not an accident            ║
  ║                                                                              ║
  ║  STEP 1  src/oinsmiles/oin/canonical_slots.py   (numpy-only header)           ║
  ║    GEOMETRY_VERTICES ─ one vertex table (closes TD-005)                       ║
  ║    derive_rotation_group(vertices)  ── filters det > 0  ◄── the enantiomer     ║
  ║    lexmin_vertex_signature(geo, vcolor) -> (signature, perm)     guard rail    ║
  ║        └─ compare.py::_polyhedron_signature takes ONLY `signature`             ║
  ║           ⇒ KEY BEHAVIOUR BYTE-IDENTICAL, proven on 13,187 capstone keys       ║
  ║             SHA-256 identical before/after                                    ║
  ║    canonical_slot_permutation / canonical_slot_map / canonicalize_oin_slots    ║
  ║                                                                              ║
  ║  STEP 2  three order-dependent SEAMS neutralized together                      ║
  ║    _align_to_pai pivot ─ np.min(candidates) ─► _canonical_pivot                ║
  ║    _align_to_pai Z-sign ─ sum(z*(i+1)^3)   ─► _canonical_z_sign                ║
  ║    get_input_order_key ─ min(atom index)    ─► donor position in canon. SMILES ║
  ║                                                                              ║
  ║  STEP 3  the relabel is a POST-PASS on the FINISHED inline string              ║
  ║    utils/xyz2mol.get_oin_string  ──►  canonicalize_oin_slots(oin)              ║
  ║    NOT inside the aligner: leaves the lex-MAX loop, the eta RC1 content        ║
  ║    swap and the heading-atom tiers (all read item["slot"]) untouched           ║
  ║                                                                              ║
  ║  OUTCOME                                                                      ║
  ║    ✔ byte-stable 88 → 100 of ~150   ✔ 0 molecules broken, either arm           ║
  ║    ✔ key-level defects 23 → 23 / 18 → 18  (the over-folding guard, at scale)   ║
  ║    ◆ acceptance target `slot_renumber → ~0` MEASURED UNREACHABLE at this seam  ║
  ║      32/32 residual pairs are `same_vcolor_identical`                          ║
  ║    ✚ ALSO SHIPPED HERE: OIN_STABLE_METAL_AC (AC valence-capping order)         ║
  ║    ‼ MEASURED 0/150 molecules emit a metal `@` tag  ⇒ LANE 5 RESCOPED          ║
  ╚═══════════════════════════════╤══════════════════════════════════════════════╝
                                  │ exports canonical_slot_permutation /
                                  │ canonical_slot_map
                                  ▼
                    ══ WAVE C: LANES 5 and 6 ══
                    (Lane 6 measured this dependency does NOT exist for P3;
                     Lane 5 was never started in v0.4.5)

Legend  ✔ delivered and measured   ✖ known-red by design   ◆ target refuted
        ✚ extra deliverable        ‼ finding that rescoped another lane
```

## Initial assumptions and hypothesis

1. **The blocking relationship is real and one-directional.** Lane 2's lex-min runs over vertex
   **colours**, and a colour is `(ligand body, donor element, winding char)`. If the body is not
   canonical the colours are not canonical, so the whole minimization is over a moving target. Lane 2
   therefore cannot start, or be measured, until Lane 1 exists — and Lane 2's only meaningful A/B arm
   is the **stacked** one.
2. **The target is `slot_renumber` 315 → ~0** (from the v0.4.4 capstone `bucket_report.json`), plus
   whatever share of `rdkit_canonical` turned out to be slot-or-order drift rather than body drift.
3. **The plan's central premise:** *"the permutation `_polyhedron_signature` achieved its lex-min
   with **is** the canonical slot relabeling."* `compare.py` had computed that lex-min colored-vertex
   signature since v0.4.4 — that is what made the comparison key fac/mer-aware — and simply **threw
   the permutation away**. Recovering it was believed sufficient.
4. **The risk is over-folding, and it is the whole risk of the lane.** A relabeling that folds past
   the geometry's own symmetry silently *merges two real isomers*, and unlike drift, a merge is
   invisible. Every guard in the lane exists for that.
5. **The charter's fragment-order key** was `(canonical body, min canonical slot, sorted canonical
   slots)`.
6. **`get_input_order_key` should be deleted.**

## What was actually found

### Confirmed

- **Key behaviour is preserved exactly.** `_polyhedron_signature` was refactored onto a shared helper
  returning `(signature, perm)`; the existing key path consumes only `signature`. Proven, not argued:
  `canonical_roundtrip_key` recomputed over **all 13,187 capstone OIN strings**, **SHA-256 identical**
  before and after. This is what made a change inside the comparison key's own machinery safe to make
  at all.
- **The rotation group was incomplete, in one specific and consequential way.**
  `oin_aligner._brute_force_symmetries` built the group by brute-forcing Euler triples from the fixed
  grid `[0, 90, 120, 180, 240, 270]`. **That grid cannot express a 72° five-fold rotation**, so on
  **PBP** (pentagonal bipyramidal) it found **2 of the 10** proper rotations — the encoder could not
  canonicalize a pentagonal-bipyramidal equatorial labeling at all. It agreed with
  `derive_rotation_group` on the other **ten** geometries and never invented a non-rotation.
  **The agreement test was written and shipped FIRST**, passing 10/11 geometries and failing PBP at
  2 vs 10, so the fix landed as a *measured improvement* rather than as an accident nobody noticed.
  `tests/unit/test_canonical_slots.py::TestAlignerSymmetriesAreUnified` pins both halves: agreement on
  every geometry, and the specific C₅ generator `(0,1,3,4,5,6,2)` present in the PBP group, so an
  angle-grid search cannot silently come back.
- **The encoder and the comparison key had been disagreeing about PBP on the default path.**
  `compare.py::_polyhedron_signature` has used the *derived* group **ungated since v0.4.4** while the
  encoder used the Euler grid. Unifying removes a live disagreement rather than introducing anything.
  Measured blast radius of the unification, encoder-vs-encoder byte diff over 150 molecules with
  levers OFF: **1/149 (0.67%)** — `IMACOO_comp_0 [Zr_PBP]` only. Geometry incidence in that sample:
  `SPL` 46, `OCT` 33, `TET` 32, `LIN` 11, `TPL` 9, `TBP` 7, `SPY` 5, `TPY` 4, **`PBP` 2**.
- **The relabeling is measurably inert at isomer level, in both directions.** Key-level defects
  23 → 23 (slots lever alone) and 18 → 18 (stacked). That is the over-folding guard at corpus scale:
  had the canonicalization merged a real isomer, the comparison key would have registered it.

### Measured — the four arms (seed 42, so every arm samples the same molecules)

`tools/canonicality_probe.py --n 150 --trials 2`, generator-free, on `swimlane/v045-lane2`. The probe
holds the molecular graph **fixed** and varies proper rotation, atom renumbering, and both, so the
correct answer is byte-**identical** — a known ground truth, not an inferred baseline.

| arm | byte-stable | drifted | `rdkit_canonical` | `slot_renumber` | `fragment_reorder` | `winding_star_drift` | `rotate` drift | key-level defects |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all levers OFF | 88/149 (59.06%) | 61 | 39 | 25 | **0** | **0** | **0** | 23 |
| `OIN_CANONICAL_SLOTS` only | 94/149 (63.09%) | 55 | 28 | 29 | **0** | **0** | **0** | 23 |
| Lane 1 levers only | 95/150 (63.33%) | 55 | 34 | 23 | **0** | **0** | **0** | 18 |
| **stacked (Lane 1 + slots)** | **100/150 (66.67%)** | 50 | **19** | 32 | **0** | **0** | **0** | **18** |

Per-molecule attribution, not just counters: OFF → slots-only is **6 fixed / 0 broken**; Lane-1-only
→ stacked is **5 fixed / 0 broken**.

*(The two OFF-side arms encode 149 rather than 150 because `OIN_STABLE_METAL_AC`, also on this
branch, recovers one input the base encoder fails on. Read the byte-stable **counts**, 88 → 100.)*

**`slot_renumber` rising while byte-stability rises is reclassification, not regression.** The probe's
taxonomy tests `slot_renumber` *before* `rdkit_canonical`, so a molecule whose body drift closes while
its slot drift remains **moves into** the finer class. Measured: 8 molecules
`rdkit_canonical → slot_renumber` on the OFF→slots step and 11 more on the Lane-1→stacked step, with
**zero** going the other way and **zero** going stable→drift.

### Refuted — the plan's central premise

> *"The permutation `_polyhedron_signature` achieved its lex-min with **is** the canonical slot
> relabeling"* is **necessary but not sufficient**, so the acceptance target `slot_renumber → ~0` is
> **not reachable at this seam.** Measured, not argued.

The `slot_renumber` residual in the stacked arm is **32 molecules**, and the mechanism is uniform:
**32/32 are `same_vcolor_identical`** — `compare._parse_vertex_colors` returns the *identical*
`{slot: colour}` map for the base string and the drifted one. That is fatal for a post-pass whose
relabeling is derived from `vcolor` alone: identical input ⇒ identical permutation ⇒ applying it to
two differently-labeled strings **preserves the difference exactly**. Not a bug in the relabeling —
a limit of its **input**.

And the reason the colours are identical is a deliberate design choice one layer down.
`_parse_vertex_colors` colours *every* donor of a ligand with that ligand's **whole** body, with no
chelate grouping, precisely so that "a swap between two same-coloured donors is invisible" (its own
docstring). That blindness is what lets true conformers collapse in the comparison key while fac/mer
stay distinct. It is **correct for a key** and **insufficient for an encoder**, which must emit one
string and therefore has to decide something the key is entitled not to care about: *which donor atom
of a chelate holds which slot integer.*

Worked example, `AGUKOD_comp_0` (Rh, SPL, a COD ligand spanning two cis vertices):

```
base: [Rh_SPL].[CH]{1}1=[CH]{1>}CC[CH]{0>}=[CH]{0}CC1.…{2}….[Cl]{3}
got : [Rh_SPL].[CH]{0}1=[CH]{0>}CC[CH]{1>}=[CH]{1}CC1.…{2}….[Cl]{3}
```

Bodies byte-identical, `vcolor` identical (`{0: (cod, C, ">"), 1: (cod, C, ">")}`), only which alkene
arm carries which integer moved. And the transposition `(0 1)` relating them is **not an element of
SPL's rotation group** (D₄ contains `(1 3)`, `(0 2)`, `(0 1)(2 3)`, … but not `(0 1)` alone), so the
group-theoretic fold cannot see it either.

**Atom-level verdict on the 32** (`tools/slot_drift_mechanism.py`, using
`CanonicalRankAtoms(breakTies=False)` symmetry classes — `breakTies=True` settles ties between
symmetry-equivalent atoms on the *input index*, i.e. the very dependence under measurement, so it is
the wrong instrument):

| verdict | n | meaning |
|---|---:|---|
| `automorphism` | **23** | the two donor atoms are in one symmetry class of their own fragment, so both labelings denote the **same molecule**; the encoder merely lacks a deterministic choice |
| `DISTINCT_donors` | **7** | slots land on **inequivalent** donors ⇒ one of the two strings is *wrong* — reported as a soundness class |
| unparsable | 2 | fragment would not sanitize (borane-cluster class), classified by nobody |

**Root cause of the 23, named:** the physical donor→vertex map comes from the 3D fit and is
order-free. What moves is **which of two symmetry-equivalent donor atoms RDKit writes first** in the
canonical SMILES — `_smilesAtomOutputOrder`, i.e. `CanonicalRankAtoms(breakTies=True)`, whose tie is
settled on the input index. The canonical *string* is invariant; the *map from string position to
atom* is not.

**The 7 `DISTINCT_donors` were then refuted by [Lane 9](LANE-09-inequivalent-donors.md)** — for all
seven the two competing strings place donors on vertices related by a **proper rotation** of the
coordination polyhedron (`identity 4 / in_rotation_group 4 / NOT_A_ROTATION 0`; `|Δrssd|` across all
8 pairs `0.0 … 1.2e-14`, exactly tied; donor count equals vertex count for all seven, so the group
test is exact rather than permissive; `key_equal = True` for 14/14 variants). Lane 9 shipped no
encoder change, because a fix would have frozen one of two equally faithful labelings.

### The finding that rescoped another lane

**0 of 150 molecules emit a metal `@` tag at all.** So Lane 5 (metal Δ/Λ) is not "un-folding a
collapse" — it is **creating a descriptor**, which is a different and larger piece of work with a
different risk profile. That measurement is why Lane 5 was rescoped and consequently **not started in
v0.4.5** (see [Wave C](WAVE-C-injectivity-descriptors.md)). Two corollaries recorded in
`docs/CANONICAL_SLOTS_v0.4.5.md` §9: the metal's own `@SP1`/`@OH1` tag *would* be a parity against a
donor ordering this post-pass renumbers, so there is nothing stale to break today — and
`normalize_oin_for_comparison` strips metal `@` stereo anyway, so no current measurement could see it
either.

## What was done

**→ [Lane 2's own report](LANE-02-canonical-slots.md); design document `docs/CANONICAL_SLOTS_v0.4.5.md`.**

### The new module: `src/oinsmiles/oin/canonical_slots.py`

Deliberately light imports — the module **header is numpy only**, no RDKit, no aligner — so both
`compare.py` (which advertises a light graph) and the encoder can import it freely. The string-level
half (from `_relabel_slots` down) needs `compare` and `inline` and imports them *inside* the
functions, keeping `compare.py → canonical_slots.py` acyclic.

| symbol | what it is |
|---|---|
| `GEOMETRY_VERTICES` | the single vertex-direction table, slot *i* == vertex *i*, for `LIN TPL SPL TET TPY TBP SPY OCT PBP` and the rest. **Closes open debt TD-005** — `utils.oin_aligner.TEMPLATE_SPECS` keeps only its `ref` vectors, which only winding needs, and `tests/unit/test_canonical_slots.py` cross-checks the two |
| `derive_rotation_group(vertices, tol=1e-3)` | enumerates the proper-rotation group by backtracking over spanning vertex sets, **filtering on `det > 0`** |
| `geometry_rotation_group(geo)` / `geometry_vertex_count(geo)` | memoized per-template accessors |
| `lexmin_vertex_signature(geo, vcolor) -> (signature, perm)` | the shared helper. `compare.py::_polyhedron_signature` now consumes **only `signature`** ⇒ key behaviour byte-identical |
| `canonical_slot_permutation(geo, vcolor) -> {int: int}` | the relabeling for one polyhedron |
| `canonical_slot_relabeling(oin_string)` / `canonical_slot_map(oin_string)` | the string-level entry points; `canonical_slot_map` is the export Lanes 5/6 must use |
| `canonicalize_oin_slots(oin_string) -> str` | the post-pass the encoder calls |

### Where the relabel is applied, and why there

**As a post-pass on the finished inline string**, in `utils/xyz2mol.get_oin_string` — **not** inside
`_permute_and_serialize`'s lex-**max** loop. This is a deliberate blast-radius decision and it is
worth understanding rather than re-litigating:

- The inline string is exactly the representation `compare._parse_vertex_colors` already reads.
  Reusing that function **verbatim** means the encoder's canonicalization and the comparison key's
  canonicalization consume the same bytes through the same code, so they agree *by construction*
  rather than by two parallel implementations staying in sync.
- It leaves three things running untouched beforehand — the geometric Kabsch fit, the **eta RC1
  content swap**, and the **heading-atom tiers** — *all* of which read `item["slot"]`. Surgery inside
  the lex-max loop would have put every one of them in the blast radius.
- **Why it is safe for winding.** Winding characters (`>` `<` `^`) are carried through **verbatim**,
  and cannot be invalidated by the relabel for two independent reasons: winding is measured against
  each ring's **own** metal→centroid axis (`oin_aligner._determine_winding`), never against the
  slot's template direction, and it is computed against the original slot's template frame **before**
  the relabel; and the group contains **only proper rotations**, which preserve circulation sense.
  That is precisely why `derive_rotation_group` filters to `det > 0`. (`_parse_vertex_colors` folds
  `^` to `>` for *colouring* only; the emitted string keeps whatever the aligner computed.)

### The tie-break, and why it is deliberately tested

Among the permutations achieving the minimal signature, `canonical_slot_permutation` takes the
**lexicographically smallest permutation tuple**. This is tested on purpose, because it is **the exact
shape of the Wave-2 trap** that destroyed the chirality it was supposed to encode: there, sorting
symmetry-equivalent axes by `ax.sign` made the axial token reflection-**invariant** (`+-` and `-+`
both rendered `-+`), and only a corpus-wide mirror audit caught it. A tie-break that consults a
*stereochemical* quantity can silently fold enantiomers; a tie-break on the permutation tuple cannot,
because the candidate set is closed under the group and the group contains only proper rotations.

There is a second, subtler consequence, and it is why **`canonical_slot_permutation` is the wrong
export for Lanes 5 and 6**: the lex-min *signature* is invariant, but the permutation achieving it
need not be **unique** — when the coloured polyhedron has a nontrivial colour-preserving rotation
stabilizer, several permutations tie, and that function breaks the tie on the permutation tuple, which
is a property of the **incoming** labeling. `canonical_slot_map` refines the tie on the **rendered
output**, which is invariant by construction (proof in its docstring: for `L' = g·L` the candidate
sets are equal, because `p ↦ p·g` is a bijection of the group). The export contract is therefore:

```python
from oinsmiles.oin.canonical_slots import canonical_slot_map
canonical_slot = canonical_slot_map(oin_string)[slot_as_it_appears_in_that_string]
```

It is the **identity** when `OIN_CANONICAL_SLOTS` is on, because the encoder has already applied it —
which is the point: a caller written against this helper is correct with the lever either way and
stays correct after promotion, and it is idempotent on an already-canonical string. What is unique is
the canonical **string**, always; the per-donor map is unique only up to the coloured polyhedron's
**rotational** automorphism group, and exactly up to it. `fac`-M(ppy)₃ has a real C₃ axis through its
three equivalent ligands, so three donor labelings are equally canonical while the emitted string
stays single-valued. That is not a defect to route around: every member of that group is a **proper**
rotation, so a descriptor that flips under reflection — metal Δ/Λ, an eta winding, an axial sign —
takes the *same* value on all of them; and a descriptor that is *not* invariant under the automorphism
is not a property of the molecule at all, only of which interchangeable ligand you chose to call
"first".

### The three order-dependent seams, neutralized together

Slot relabeling alone is not enough, because two upstream steps make the *frame itself* a function of
the input numbering, and a different frame produces a different geometric fit:

| seam | was | now (lever ON) |
|---|---|---|
| `_align_to_pai` pivot | `np.min(candidates)` — lowest **file** index among atoms tied for max distance from the metal | `_canonical_pivot`: `(-mass, sorted multiset of interatomic distances)` — an atomic property plus a rigid-motion invariant |
| `_align_to_pai` Z-sign | `sum(z_i * (i+1)**3)` — a weight built from the atom's **position in the file** | `_canonical_z_sign`: mass-weighted **odd** moments in z, tried in order; plain sums, so unchanged by renumbering, and odd, so they actually decide `+z` vs `−z` |
| `get_input_order_key` | `min(original atom indices)` | donor position inside the fragment's own canonical SMILES (`s_idx`, off `_smilesAtomOutputOrder`) |

### Fragment order, and the charter deviation

Re-derived from the canonical slots as `(minimum canonical slot, fragment text)`, with the metal
fragment pinned first — `fragments[0]` is a load-bearing project invariant that the generator, the
inline parser and the comparison key all assume. `minimum canonical slot` is by itself a **total**
order on the coordinated fragments, because a slot belongs to exactly one fragment; the text term only
ever separates two *uncoordinated* fragments.

The charter proposed `(canonical body, min canonical slot, sorted canonical slots)`. Body-primary was
**rejected on purpose**: it is equally invariant but not *more* invariant, and it reorders fragments
relative to the pre-existing convention (step 7 of `get_oin_string` already sorted on `_sort_slot`
primary) for no gain — which would multiply byte churn at promotion time across every molecule rather
than only the ones whose slots actually move. `sorted canonical slots` is unreachable: two fragments
cannot share a minimum slot.

### `get_input_order_key` is superseded, not deleted

Deleting it outright would change levers-OFF bytes, which the release's default-OFF discipline
forbids — and the lane's brief asked for both, a contradiction it resolved and documented
(`docs/CANONICAL_SLOTS_v0.4.5.md` §6) rather than silently picking one. The supersession is total on
the canonical path: its **second** use (step 7, the final fragment sort) is overwritten wholesale by
the post-pass; its **first** use (step 5, before the aligner runs) only ever separates fragments that
already tied on mass, binder mass *and* body SMILES — i.e. two copies of one ligand, which get
identical vertex colours, so swapping them leaves `vcolor` unchanged and the lex-min labeling
identical. The empirical form of the same claim is the acceptance test
`tests/unit/test_canonical_slot_invariance.py`: permute the input atom order, the emitted string must
be byte-identical.

### Also shipped on this branch: `OIN_STABLE_METAL_AC`

Not in the lane's charter. `swimlane/v045-lane2` @ `8bf9df61` closes `DUDREA_comp_0`, whose geometry
tag flipped `[Y_SPY]` → `[Y_TET]` under renumbering. The only order-dependent step in adjacency
perception is the **valence-capping loop**, which iterated in atom-index order: capping the metal
before vs after a bridging hydride decided whether the Y–H bond survived (degree 5/`SPY` vs
4/`TET`). Measured: the adjacency matrix differed in **3 of 8** renumberings with the lever off,
**0 of 8** with it on. Note where this lands — the `structural` bucket, 16.1% of the gap to 100%,
which the handoff had written off as out-of-scope perception research. `DUDREA_comp_0` was never a
slot problem; do not reopen it as one.

### The goldens that moved, and the proof each is a relabeling

With `OIN_CANONICAL_SLOTS` on, three of the six pinned fixtures move and three are already canonical.
**Every moved string has an unchanged `canonical_roundtrip_key`** — that is the assertion which would
catch a canonicalization that merged two isomers, and it was checked *before* re-pinning, not after.

| fixture | off | on | key |
|---|---|---|---|
| CisPlatin | `[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}` | `[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}` | identical |
| TransPlatin | `[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}` | `[Pt_SPL].N{0}.[Cl]{1}.N{2}.[Cl]{3}` | identical |
| fac-Ir(ppy)₃ | `…c{0}…n{3}1.c{5}…n{1}1.c{2}…n{4}1` | `…c{0}…n{5}1.c{2}…n{1}1.c{4}…n{3}1` | identical |
| PdCl₂-R-BINAP | — | P moves to slots 2,3 | identical |
| mer-Ir(ppy)₃, Ferrocene, Cis-PtCl₂(en) | — | unchanged | — |

**Why they moved, in one sentence:** the lex-min runs over vertex **colours** and `"N" < "[Cl]"`
bytewise, so the amines take the low slots. **Cis is still cis** — the chlorides land on 2 and 3,
which are adjacent. All of it is pinned by `TestLeverOnGoldens`, which asserts the moved strings *and*
key-invariance *and* the unmoved set, so gratuitous churn is caught too; that turns the Wave D
promotion A/B into a diff against a committed expectation rather than a judgement call.

### The wider fold, specified and declined

`docs/CANONICAL_SLOTS_v0.4.5.md` §7a specifies the fix that would close the 23 `automorphism`
residuals: fold, in addition to the rotation group, the permutations exchanging donors **within one
fragment** that are (a) in the same `breakTies=False` symmetry class and (b) the same colour. Scoped
that narrowly it does not touch metal Δ/Λ (which lives in the arrangement *across* fragments, still
folded only by proper rotations) and it cannot touch inequivalent donors (different classes → never
folded).

**It was declined, twice, and the second decision is the durable one.** Lane 2 gave three reasons:
it folds past the geometry's own symmetry, it needs a Δ/Λ guard that does not exist, and it would
over-fold the 7 `DISTINCT_donors`. Lane 9 **refuted the third** — all seven are already related by a
proper rotation, which the post-pass folds over anyway. The release still declined, on stronger
grounds:

- **the payoff is zero at the headline** — `slot_renumber` drift lands in `key_equal`, which *already
  counts as passing*, so widening the fold cannot move the round-trip number by a single molecule;
- **the risk is the worst failure mode in the release** — over-folding silently merges two real
  isomers, and unlike drift a merge is invisible;
- **Lane 9 identified a better fix.** The real defect behind all 30 residuals is one sentence: *the
  choice among group-related labelings is settled on atom order, in three separate places* — the
  Kabsch tie (`_map_to_template`'s strict `<`), `_permute_and_serialize`'s lex-max, and the body's
  atom output order. Making those three invariant is narrower, safer and principled. Note a fix aimed
  at the Kabsch tie **alone closes at most 3 of 7** — for 4 of them the fit never moved (π = identity).

## Dead ends and refutations

| tried / believed | what killed it |
|---|---|
| "the argmin permutation of `_polyhedron_signature` **is** the canonical relabeling" (the plan's premise) | necessary, not sufficient. **32/32** residual pairs are `same_vcolor_identical` — identical `vcolor` input ⇒ identical permutation ⇒ the difference survives the post-pass exactly |
| acceptance target `slot_renumber → ~0` | unreachable at this seam. `slot_renumber` in fact **rises** (25 → 32) while byte-stability rises — reclassification, verified per-molecule as 19 moves `rdkit_canonical → slot_renumber`, 0 the other way, 0 stable→drift |
| apply the relabel inside `_permute_and_serialize`'s lex-max loop | rejected on blast radius: the lex-max loop, the eta RC1 content swap and the heading-atom tiers all read `item["slot"]`. A post-pass on the finished string reuses `compare._parse_vertex_colors` verbatim instead, so encoder and key agree by construction |
| the charter's `(canonical body, min canonical slot, sorted canonical slots)` fragment key | body-primary is not *more* invariant (min-slot is already total on coordinated fragments) and multiplies promotion-time byte churn; `sorted canonical slots` is unreachable since two fragments cannot share a minimum slot |
| delete `get_input_order_key` | would break levers-OFF byte-identity, which acceptance also required. Superseded instead, and the contradiction in the brief documented rather than silently resolved |
| gate the rotation-group unification behind the lever | it would preserve a **live encoder/key disagreement** — `compare.py` has used the derived group ungated since v0.4.4. Left unconditional, with its 1/149 blast radius measured and flagged rather than decided silently |
| `CanonicalRankAtoms(breakTies=True)` as the symmetry-class instrument for the residual triage | it settles ties between symmetry-equivalent atoms on the *input index* — the very dependence under measurement. `breakTies=False` used instead (with `includeChirality` defaulting True, so two constitutionally-equivalent branches with different configurations correctly land in different classes) |
| "**7 molecules have slots on inequivalent donors, so one string is simply wrong**" — reported and acted on | **refuted by [Lane 9](LANE-09-inequivalent-donors.md).** All seven pairs are related by a proper rotation of the polyhedron; `|Δrssd|` 0.0 … 1.2e-14; `key_equal` 14/14. The classifier asked a per-fragment, positionally-zipped question about a **global** relabeling (`RUBTIS`'s COD arms swapped too, making the relabeling `(0 1)(2 3) ∈ D₄`), plus two mechanical bugs (`zip` assumed a stable fragment order that the post-pass derives *from* slots; haptic donors keyed on the marker's first occurrence). Renamed `distinct_donors_LOCAL` with a warning |
| take §7a's wider fold now that Lane 9 removed one objection | declined: zero headline payoff, worst-in-release risk profile, and a narrower alternative exists |
| assume the Kabsch tie-break is where the residual drift enters | for **4 of the 7** the fit did not move at all (π = identity); the drift enters downstream in `_permute_and_serialize`'s lex-max, its homogeneous-sort rank assignment, and which of two equivalent donors the body renders first |

**One finding in passing, unowned:** `IMACOO`'s tied set is 2 of PBP's 10 group elements with the
nearest non-tied candidate only **6.8e-08** away — the stored 7-decimal C₅ coordinates split the
orbit, so **PBP ties break on rounding noise.** A real, separate order-sensitivity. And 2 of the 32
residuals are `unparsable` (borane class), classified by nobody.

## Where it landed

- **Branch** `swimlane/v045-lane2`, tip **`840eab84`**, 6 commits:
  `cc0758d7` (extract `oin/canonical_slots`, return the lex-min permutation) ·
  `51770c93` (merge `swimlane/v045-lane1` in as the shared root — the hard dependency, made explicit
  in the history) · `3722b18e` (`WIP(lane2)` — the mid-task hard stop) ·
  `8bf9df61` (`OIN_STABLE_METAL_AC`) · `12898c5a` (measure the lever; the plan's premise is
  insufficient) · `840eab84` (pin exactly which goldens move, and that each is key-identical).
- **Merged into `release/v0.4.5`** as **`c63d3404`**, third in the land order (after Lane 7 and
  Lane 1), then into local `main` as part of **`0d165845`**, tag **`v0.4.5`**. Not pushed.
- **Levers:** `OIN_CANONICAL_SLOTS` and `OIN_STABLE_METAL_AC` — both **default-ON** since v0.4.5 via
  `src/oinsmiles/oin/levers.py::_DEFAULT_ON`, on the evidence in `docs/PROMOTION_GATE_v0.4.5.md`.
- **Code:** `src/oinsmiles/oin/canonical_slots.py` (new) ·
  `src/oinsmiles/oin/compare.py::_polyhedron_signature` (refactored onto the shared helper, key
  behaviour byte-identical) · `src/oinsmiles/utils/oin_aligner.py`
  (`_brute_force_symmetries` now delegates; `_canonical_pivot`, `_canonical_z_sign`) ·
  `src/oinsmiles/utils/xyz2mol.py::get_oin_string` (the post-pass) ·
  `src/oinsmiles/utils/xyz2mol_local.py` (the capping-order fix).
- **Guards:** `tests/unit/test_canonical_slots.py` (20 new tests, incl.
  `TestAlignerSymmetriesAreUnified` and `TestLeverOnGoldens`) ·
  `tests/unit/test_canonical_slot_invariance.py` (permute the input, demand byte-identity) ·
  `tests/unit/test_facmer_key.py` and `tests/unit/test_isomer_divergence.py` (the over-folding
  vetoes — those two guards **are** the entire risk of this lane) ·
  `tests/unit/test_regression_stability.py` (the default-path golden pin).
- **Suite on the branch:** **674 OK / 3 skip / 3 xfail**, from its own 668 baseline; ruff clean.
- **Debt closed:** **TD-005** — one home for the geometry vertex table
  (`oin/canonical_slots.py::GEOMETRY_VERTICES`), replacing the duplicated `TEMPLATES`/vertex data.
- **Proof of non-regression on the key:** `canonical_roundtrip_key` over all **13,187** capstone OIN
  strings, **SHA-256 identical** before and after the `_polyhedron_signature` refactor.
- **A known promotion-time consequence, flagged deliberately rather than papered over:**
  `tests/unit/test_regression_stability.py` pins the **default path** and does not clear the
  environment, so running the suite with `OIN_CANONICAL_SLOTS` exported fails exactly three goldens
  (`test_cisplatin`, `test_transplatin`, `test_fac_irppy3`) and nothing else. Left as-is by the lane
  because that file is the project's shared default-path pin; resolved in
  [Wave D](WAVE-D-integrate-promote-release.md)'s triage.
- **Design document:** `docs/CANONICAL_SLOTS_v0.4.5.md`. Lane report:
  [LANE-02-canonical-slots.md](LANE-02-canonical-slots.md).

## Open questions / for the next agent

1. **The 23 `automorphism` residuals are still open, and the recommended fix is *not* §7a's wider
   fold.** Make invariant the three places where the choice among group-related labelings is settled
   on atom order: `_map_to_template`'s strict `<` (the Kabsch tie), `_permute_and_serialize`'s lex-max,
   and the body's atom output order. A fix aimed only at the Kabsch tie closes **at most 3 of 7** of
   the Lane-9 set, so all three sites are needed. If §7a is ever revisited, it needs a Δ/Λ
   tris-chelate divergence fixture as a minimum guard — `tests/fixtures/ZUMNEC.xyz` now exists for
   exactly that (see [Lane 7](LANE-07-research-residuals.md)).
2. **`_parse_vertex_colors` is blind by design, and that blindness is now the binding constraint.**
   Any attempt to close `slot_renumber` further has to give the encoder a *finer* colouring than the
   key's, without giving the key one — because the key's blindness is what lets true conformers
   collapse while fac/mer stay distinct. Do not "fix" `_parse_vertex_colors` in place.
   Per-eta-ring colour (`_eta_automorphism_class`) is the one refinement already specified, and it is
   [Lane 3](LANE-03-winding-residual.md)'s ask, for QIGZAJ.
3. **PBP ties break on rounding noise.** `IMACOO`'s nearest non-tied candidate is 6.8e-08 away
   because the stored C₅ coordinates carry 7 decimals. Either raise the template precision or make
   the tie-break exact; measure on the PBP population (2 of 150 in the probe sample) rather than on
   `IMACOO` alone.
4. **2 of the 32 residuals are `unparsable` (borane class) and no lane owns them.** With
   `OIN_BORON_CAGE` now default-ON in v0.4.6 they may have become parsable — re-run
   `tools/slot_drift_mechanism.py` and re-triage.
5. **`canonical_slot_permutation` vs `canonical_slot_map` is a live footgun.** Any future descriptor
   lane must use `canonical_slot_map(oin_string)`; `canonical_slot_permutation(geo, vcolor)` breaks
   its tie on the **incoming** labeling and is not invariant when the coloured polyhedron has a
   nontrivial colour-preserving stabilizer. This is documented in the module and in
   `docs/CANONICAL_SLOTS_v0.4.5.md` §8, and it has already been mis-cited once.
6. **The metal `@` tag is now reachable.** `0/150` was true at Lane 2's measurement because no
   descriptor existed; [Lane 5](LANE-05-metal-delta-lambda-P1.md) built one in v0.4.6 behind
   `OIN_EMIT_METAL_CONFIG`. If that lever is ever promoted, re-check §9's claim that there is "nothing
   stale to break" — a metal parity tag *is* a parity against a donor ordering this post-pass
   renumbers, and `normalize_oin_for_comparison` strips metal `@` stereo, so no existing measurement
   would see a breakage.
