# Lane 2 — Canonical coordination slots

**What this lane was for:** make each donor's `{n}` slot integer a **graph invariant** instead of
"whichever idealized template vertex its 3D direction vector happened to fit best, with ties
settled on the raw input atom index" — and re-derive the emitted fragment order from those
canonical slots.

---

## ELI5

An OIN string labels each atom that touches the metal with a small integer — its "slot" — meaning
*which corner of the coordination polyhedron* it sits on. Those integers used to come from fitting
the measured 3D directions onto an idealized shape, and when two corners were equally good the
code picked whichever atom appeared earlier in the input file. So re-ordering the lines of an XYZ
file (same molecule, same shape) could hand out different slot numbers and produce a different
string. The fix: keep the *shape* from the 3D fit — that is what says cis vs trans and fac vs mer —
but then, among only the corner-relabelings the shape's own symmetry allows, pick the alphabetically
smallest one. "Only what the symmetry allows" is the whole trick: sorting ligands by priority and
stamping them onto 0,1,2,… would erase cis/trans, because cis/trans *is* which corners equal
ligands occupy. And only **rotations** are allowed, never mirror images, because folding in a
mirror would merge a molecule with its non-superimposable twin.

## The work, visually

```
                       3D fit to idealized template (UNCHANGED — decides the ISOMER)
                                        │
                                        ▼
        ┌─────────── frame-level order dependence, neutralized together ───────────┐
        │  xyz2mol.py::_align_to_pai                                               │
        │    pivot     : np.min(candidates)  ── FILE index ──┐                      │
        │                → _canonical_pivot: (-mass, sorted multiset of            │
        │                  interatomic distances)             │ atomic property   │
        │    z-sign    : sum(z_i * (i+1)**3) ── FILE position ─┤ + rigid-motion    │
        │                → _canonical_z_sign: mass-weighted    │ invariant         │
        │                  ODD moments in z (plain sums)      ─┘                   │
        │  xyz2mol.py::get_input_order_key                                         │
        │    min(original atom indices) → donor position inside the fragment's own │
        │    canonical SMILES (s_idx, off _smilesAtomOutputOrder)                  │
        └──────────────────────────────────────────────────────────────────────────┘
                                        │
                       inline OIN string built as before (slots, winding, eta)
                                        │
        ┌──────── POST-PASS: canonical_slots.canonicalize_oin_slots (OIN_CANONICAL_SLOTS) ────────┐
        │                                                                                        │
        │  compare._parse_vertex_colors(normalize_oin_for_comparison(body))   ← VERBATIM reuse    │
        │        └─► vcolor = {slot: (fragment_body, donor_element, winding)}                     │
        │                                                                                        │
        │  group = geometry_rotation_group(geo)      ← derive_rotation_group, det > 0 only        │
        │                                              OCT 24, PBP 10, unknown geo → identity    │
        │                                                                                        │
        │  for perm in group:            (group arrives SORTED ⇒ deterministic argmin)            │
        │      arr        = colours placed at perm[slot]                                          │
        │      out, ord   = _render(frags, metal_pos, mapping)   ← metal pinned fragments[0]      │
        │      key        = (tuple(arr),      ord,               perm)                            │
        │                    ^^^^^^^^^^       ^^^                ^^^^                              │
        │            1. lex-min signature   2. rendered      3. tie-break on the                  │
        │               = compare.py's        output           permutation tuple only              │
        │               _polyhedron_signature (invariant:      (never on colour content —          │
        │               ⇒ encoder & key         p↦p·g is a      a Y2-style sign sort made a         │
        │               stay in lockstep        bijection)      token reflection-invariant)        │
        │      keep min(key)                                                                      │
        │                                                                                        │
        │  winding chars > < ^ carried through VERBATIM (measured against each ring's own          │
        │  metal→centroid axis, and proper rotations preserve circulation sense)                  │
        │  |ax:±| suffix split off by _AXIAL_SUFFIX_RE and re-appended unchanged                  │
        └────────────────────────────────────────────────────────────────────────────────────────┘
                                        │
                            fragment order = (min canonical slot, fragment text)

    UNCONDITIONAL (not behind the lever):
      oin_aligner._brute_force_symmetries  ──delegates to──►  canonical_slots.derive_rotation_group
      Euler grid [0,90,120,180,240,270] found PBP 2/10 rotations (no 72° 5-fold) → now 10/10
      blast radius measured: 1 of 149 molecules changes bytes (IMACOO_comp_0 [Zr_PBP])
```

Legend: everything in the POST-PASS box and the frame-level box is behind `OIN_CANONICAL_SLOTS`
(default-**ON** since v0.4.5); the rotation-group unification is deliberately unconditional.

## Initial assumptions and hypothesis

1. **The plan's premise:** "the permutation that `compare._polyhedron_signature` achieved its
   lex-min with **is** the canonical slot relabeling." `compare.py` had computed that lex-min
   colored-vertex signature since v0.4.4 (it is what makes the round-trip key fac/mer-aware) and
   simply *discarded* the achieving permutation. Recovering it was expected to be most of the work.
2. **The acceptance target** was `slot_renumber` → ~0 (315 molecules on the v0.4.4 capstone sweep).
3. **The charter proposed** a fragment order of `(canonical body, min canonical slot, sorted
   canonical slots)`.
4. **Belief going in:** slot relabeling alone would be enough; the geometric fit itself was assumed
   presentation-independent.

## What was actually found

**Confirmed — the refactor is behaviour-preserving.** `_polyhedron_signature` was re-based onto a
shared helper `canonical_slots.lexmin_vertex_signature` returning `(signature, perm)`; the key path
takes only `signature` and drops `perm`. Proof: `canonical_roundtrip_key` recomputed over all
**13,187 capstone OIN strings**, **SHA-256 identical before and after**
(`docs/V045_STATUS_2026-07-25.md`).

**Confirmed — the encoder and the key disagreed about the rotation group.**
`oin_aligner._brute_force_symmetries` brute-forced Euler triples from the fixed grid
`[0, 90, 120, 180, 240, 270]`, which **cannot express a 72° five-fold rotation**. On PBP it found
**2 of the 10** proper rotations, so the encoder could not canonicalize a pentagonal-bipyramidal
equatorial labeling at all — while `compare.py` had used the *derived* group, ungated, since
v0.4.4. It agreed on the other ten geometries and never invented a non-rotation, which is what
made unifying safe. Measured blast radius of unifying, deterministic encoder-vs-encoder byte diff
over 150 molecules with levers OFF: geometry incidence `SPL` 46, `OCT` 33, `TET` 32, `LIN` 11,
`TPL` 9, `TBP` 7, `SPY` 5, `TPY` 4, **`PBP` 2**; **byte-different 1/149 (0.67%)** —
`IMACOO_comp_0 [Zr_PBP]` only.

**Refuted (partly) — slot relabeling alone is not enough.** Two upstream seams make the *frame*
itself a function of the input numbering, and a different frame gives a different geometric fit:
`_align_to_pai`'s pivot (`np.min(candidates)` = lowest file index among atoms tied for max distance
from the metal) and its Z-sign metric (`sum(z_i * (i+1)**3)`, a weight built from the atom's
position in the file). Live evidence: `DUDREA_comp_0` flips its geometry classification
**`[Y_SPY]` → `[Y_TET]` under pure renumbering**.

**Measured — `tools/canonicality_probe.py --n 150 --trials 2`, seed 42, branch
`swimlane/v045-lane2`.** Lane 2 declares a hard dependency on Lane 1 (the canonical ligand body
*is* the colour the slot lex-min minimizes over), so the stacked arm is the one that matters:

| arm | byte-stable | drifted | `rdkit_canonical` | `slot_renumber` | `fragment_reorder` | `winding_star_drift` | `rotate` drift | key defects |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all levers OFF | 88/149 (59.06%) | 61 | 39 | 25 | **0** | **0** | **0** | 23 |
| `OIN_CANONICAL_SLOTS` only | 94/149 (63.09%) | 55 | 28 | 29 | **0** | **0** | **0** | 23 |
| Lane 1 levers only | 95/150 (63.33%) | 55 | 34 | 23 | **0** | **0** | **0** | 18 |
| **stacked (Lane 1 + slots)** | **100/150 (66.67%)** | 50 | **19** | 32 | **0** | **0** | **0** | **18** |

(The OFF-side arms encode 149 not 150 because `OIN_STABLE_METAL_AC` — this branch's own AC fix —
also recovers one input the base encoder fails on. Read the counts, 88 → 100.)

Per-molecule attribution: **OFF → slots only = 6 fixed, 0 broken**; **Lane 1 → stacked = 5 fixed,
0 broken**. `slot_renumber` *rising* is **reclassification, not regression** — the probe's taxonomy
tests `slot_renumber` before `rdkit_canonical`, so a molecule whose body drift closes while its
slot drift remains *moves into* the finer class: 8 molecules `rdkit_canonical → slot_renumber` on
the OFF→slots step, 11 more on the Lane1→stacked step, **zero** the other way, **zero**
stable→drift. **The slot lever moves nothing at isomer level in either direction** (key defects
23 → 23 alone, 18 → 18 stacked) — the corpus-scale over-folding guard.

**Refuted — the acceptance target is unreachable at this seam.** `slot_renumber` residual in the
stacked arm is **32 molecules**, and the mechanism is uniform: **32/32 are
`same_vcolor_identical`**, i.e. `compare._parse_vertex_colors` returns the *identical*
`{slot: colour}` map for the base string and the drifted one. A post-pass that derives its
relabeling from `vcolor` alone therefore computes the *same* permutation for both strings and
**preserves the difference exactly**. Not a bug in the relabeling — a limit of its input.
`_parse_vertex_colors` colours *every* donor of a ligand with that ligand's **whole** body, with no
chelate grouping, deliberately, so "a swap between two same-coloured donors is invisible" (its own
docstring). Correct for a key; insufficient for an encoder, which must decide *which donor atom of
a chelate holds which slot integer*.

Worked example, `AGUKOD_comp_0` (Rh, SPL, a COD ligand spanning two cis vertices):

```
base: [Rh_SPL].[CH]{1}1=[CH]{1>}CC[CH]{0>}=[CH]{0}CC1.…{2}….[Cl]{3}
got : [Rh_SPL].[CH]{0}1=[CH]{0>}CC[CH]{1>}=[CH]{1}CC1.…{2}….[Cl]{3}
```

Bodies byte-identical; `vcolor` identical (`{0: (cod, C, ">"), 1: (cod, C, ">")}`); only which
alkene arm carries which integer moved. And the transposition `(0 1)` relating them is **not an
element of SPL's rotation group** (D4 contains `(1 3)`, `(0 2)`, `(0 1)(2 3)`, … but not `(0 1)`
alone), so widening the group-theoretic fold cannot reach it either.

**Atom-level verdict on the 32** (`tools/slot_drift_mechanism.py`, using
`CanonicalRankAtoms(breakTies=False)` symmetry classes — `breakTies=True` settles ties on the input
index, the very dependence under measurement, so it is the wrong instrument; `includeChirality`
defaults True, so two constitutionally-equivalent branches with different configurations correctly
land in different classes):

| verdict | n | meaning |
|---|---:|---|
| `automorphism` | **23** | both donor atoms in one symmetry class of their own fragment ⇒ both labelings denote the **same molecule**; the encoder merely lacks a deterministic choice |
| `DISTINCT_donors` | **7** | slots land on **inequivalent** donors ⇒ one of the two strings is **wrong**: a soundness defect, not a canonicality one. `IMACOO_comp_0`, `RUBTIS_comp_0`, `VEXHIR_comp_0`, `ZACFER_comp_0`, `HAVGIW_comp_0`, `KISQAG_comp_0`, `ZOSNUS_comp_0` |
| unparsable | 2 | fragment would not sanitize (borane-cluster class) |

Root cause of the 23, named: the physical donor→vertex map comes from the 3D fit and is order-free.
What moves is **which of two symmetry-equivalent donor atoms RDKit writes first in the canonical
SMILES** — `_smilesAtomOutputOrder`, i.e. `CanonicalRankAtoms(breakTies=True)`, whose tie is
settled on the input index. The canonical *string* is invariant; the *map from string position to
atom* is not.

**Measured — nothing to break at the metal today: 0 of 150 molecules emit a metal `@` tag at all.**
The metal's own `@SP1`/`@OH1` would be a parity against a donor ordering this post-pass renumbers,
so this was checked before shipping. It is also why Lane 5 was scoped as **creating** a metal
Δ/Λ descriptor rather than un-folding an existing one. (`normalize_oin_for_comparison` strips metal
`@` anyway, so no current measurement could see it either.)

**Measured — which shipped goldens move, and that each is a relabeling.** With
`OIN_CANONICAL_SLOTS` alone, three of the six pinned `test_regression_stability.py` fixtures move
and three do not. **Every moved string has an unchanged `canonical_roundtrip_key`:**

| fixture | off | on |
|---|---|---|
| CisPlatin | `[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}` | `[Pt_SPL].N{0}.N{1}.[Cl]{2}.[Cl]{3}` |
| TransPlatin | `[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}` | `[Pt_SPL].N{0}.[Cl]{1}.N{2}.[Cl]{3}` |
| fac-Ir(ppy)₃ | `c{0}…n{3}1.c{5}…n{1}1.c{2}…n{4}1` | `c{0}…n{5}1.c{2}…n{1}1.c{4}…n{3}1` |
| mer-Ir(ppy)₃, Ferrocene, Cis-PtCl₂(en) | — | unchanged |

Cisplatin is the clearest case: the old labeling put the chlorides on 0,1 only because the fragment
sort ranks by descending mass. The lex-min runs over vertex **colours** and `"N" < "[Cl]"`
bytewise, so the amines take the low vertices. **Cis is still cis** — the chlorides land on 2 and
3, which are adjacent.

`PdCl2-R-BINAP` also moved, to
`[Pd_SPL].[Cl]{0}.[Cl]{1}.c1ccc(P{2}(…)…P{3}(…)…)cc1` (P donors on slots 2,3), likewise verified
`canonical_roundtrip_key`-identical. ⚠ **Discrepancy worth knowing:** `docs/CANONICAL_SLOTS_v0.4.5.md`
§7 says "three of the six pinned fixtures move" — that count is over
`test_regression_stability.py`'s six goldens, which do **not** include BINAP. BINAP's re-pin lives
in `tests/unit/test_stereo_roundtrip_diagnostics.py::TestEtaRingCanonicalization::test_non_eta_fragment_order_is_inert`,
whose docstring records five goldens moving at promotion (four by `OIN_CANONICAL_SLOTS`, plus
BDPP/BDNN by `OIN_STABLE_STEREO`). Both statements are true of their own fixture set; neither is
the whole picture.

## What was done

**New module `src/oinsmiles/oin/canonical_slots.py`.**

- `GEOMETRY_VERTICES` — the single vertex-direction table for 11 geometry tags (LIN, TPL, SPL, TET,
  TPY, TBP, SPY, OCT, PBP, SQA, TCT), mirroring `oin_aligner.TEMPLATE_SPECS`' `pos` entries. This
  closes open debt **TD-005** (the table was duplicated); `TEMPLATE_SPECS` keeps its `ref` vectors,
  which only winding needs, and `tests/unit/test_canonical_slots.py::TestVertexTableIsSingleSourceOfTruth`
  cross-checks the two.
- `derive_rotation_group(vertices, tol=1e-3)` — permutations preserving the full Gram
  (pairwise-dot) matrix, backtracking with a rotation-invariant per-vertex fingerprint to prune.
  For a rank-3 (spanning) vertex set the realizing orthogonal map is unique, so a permutation is
  kept iff `det > 0`; for planar/linear sets (SPL, TPL, LIN) the out-of-plane direction is free, so
  every Gram-preserving permutation extends to a proper 3D rotation and all are kept. Returns a
  **sorted** list — sortedness is what makes the lex-min argmin deterministic.
- `geometry_rotation_group(geo)` (memoized), `geometry_vertex_count(geo)`,
  `lexmin_vertex_signature(geo, vcolor) -> (signature, perm)`,
  `canonical_slot_permutation(geo, vcolor)`, `canonical_slot_relabeling(oin_string)`,
  `canonical_slot_map(oin_string)`, `canonicalize_oin_slots(oin_string)`. `VERTEX_SENTINEL =
  ("~","~","~")` sorts above every ASCII letter and bracket, so an occupied vertex always beats an
  empty one.
- Import graph deliberately light: the module header is **numpy only**; `compare` and `inline` are
  imported *inside* the string-level functions so `compare.py → canonical_slots.py` stays acyclic.

**`compare.py` re-based, not copied.** `compare._polyhedron_signature` is now a thin wrapper over
`lexmin_vertex_signature` that drops the permutation — one derivation, two consumers.

**Where the post-pass is applied, and why there.** `utils/xyz2mol.get_oin_string` (~line 2033), as
a post-pass on the **finished inline string**, *not* inside `_permute_and_serialize`'s lex-max loop.
The inline string is exactly the representation `compare._parse_vertex_colors` already reads, so
reusing that function **verbatim** makes the encoder's canonicalization and the key's agree *by
construction* rather than by two implementations staying in sync. It also leaves the geometric fit,
the eta RC1 content swap and the heading-atom tiers — all of which read `item["slot"]` — running
untouched beforehand.

**Design choices, and the alternatives rejected.**

- **Proper rotations only.** An improper operation maps a structure to its mirror image, so folding
  over reflections would collapse enantiomers and flip the eta winding sense.
  `derive_rotation_group` filters spanning sets on `det > 0` for exactly that reason. This is the
  guard rail the Y2 axial wave lacked when a tie-break sorted on a stereochemical *sign* and
  silently made the token reflection-invariant.
- **Minimize over the rendered output, not the signature alone.** Several permutations can achieve
  the minimal signature — precisely when the coloured polyhedron has a nontrivial colour-preserving
  rotation stabilizer — and picking between them on the *incoming* labeling is the input dependence
  the lane exists to remove. The sort key is `(signature, rendered_order, perm)`: the signature
  stays **primary** so the emitted labeling achieves `_polyhedron_signature`'s lex-min (re-running
  the key on the emitted string finds the identity already optimal), the rendered term is invariant
  by construction (for `L' = g·L` the candidate sets are equal because `p ↦ p·g` is a bijection of
  the group), and the third term only ever separates permutations that already render identically.
- **`canonical_slot_permutation` is explicitly NOT the function to ask "what is donor *d*'s
  canonical slot?"** — its tie-break is on the permutation tuple, a property of the incoming
  labeling. Lanes 5/6 must use `canonical_slot_map(oin_string)[slot]`, which refines the tie on the
  rendered output and is the **identity** when the lever is on, so a caller is correct either way
  and stays correct after promotion.
- **Fragment order = `(min canonical slot, fragment text)`**, metal pinned first. This is a
  deliberate deviation from the charter's `(canonical body, min canonical slot, sorted canonical
  slots)`: body-primary is equally invariant but not *more* invariant (min-slot is already a
  **total** order on coordinated fragments, since a slot belongs to exactly one fragment), and it
  would reorder fragments relative to the pre-existing convention (`_sort_slot` primary) for no
  gain, multiplying byte churn at promotion time. `sorted canonical slots` is unreachable — two
  fragments cannot share a minimum slot.
- **`get_input_order_key` is superseded, not deleted.** Deleting it outright would change
  levers-OFF bytes, which the release's default-OFF discipline forbade. Its step-7 use is
  overwritten wholesale by the post-pass; its step-5 use only ever separates fragments already tied
  on mass, binder mass **and** body SMILES (two copies of one ligand), which get identical vertex
  *colours*, so swapping them leaves `vcolor` unchanged. Under the lever its value is a tuple of
  donor positions inside the fragment's own canonical SMILES, not a file index at all.
- **The rotation-group unification is left unconditional, deliberately.** `compare.py` already
  enumerated all 10 PBP rotations while the encoder enumerated 2 — a live encoder/key
  disagreement. Unifying removes it rather than introducing anything; it folds only proper
  rotations, so no isomer information can be lost, and every pinned golden stays byte-identical.
  Flagged in the doc as the one line to gate if strict corpus-wide levers-OFF byte-identity is ever
  required.
- **Unknown geometry tags fold with the identity only** — conservative by design: an unknown
  geometry can *over*-split (miss a benign rotation) but can never wrongly merge two isomers, so it
  cannot reintroduce fac/mer blindness.

**The three frame seams**, all behind the lever: `xyz2mol._canonical_pivot` (key
`(-mass, tuple(sorted(rounded interatomic distances)))` — an atomic property plus a rigid-motion
invariant, distances rounded to 1e-6 Å so float noise from the rotation cannot reorder coincident
shells); `xyz2mol._canonical_z_sign` (mass-weighted **odd** moments in z, tried in order:
`Σm·z`, `Σm·z³`, `Σm·z·r_xy²`, `Σz`, `Σz³` — plain sums, so renumbering-invariant, and odd, so they
actually decide between `+z` and `−z`; falling through all five means every odd moment vanishes,
i.e. the structure is mirror-symmetric across the xy-plane and either orientation is correct);
and `get_input_order_key` as above.

## Dead ends and refutations

| tried | what killed it |
|---|---|
| "the lex-min permutation **is** the canonical slot relabeling" (the plan's premise) | necessary but **not sufficient**. 32/32 residual `slot_renumber` pairs are `same_vcolor_identical`, so the same permutation is computed for both strings and the difference survives. Acceptance target `slot_renumber → ~0` is **not reachable at this seam** |
| canonicalize slot labels only, leaving the frame alone | `DUDREA_comp_0` flips `[Y_SPY]` → `[Y_TET]` under pure renumbering — a different frame is a different fit. Required `_canonical_pivot` + `_canonical_z_sign` |
| charter fragment order `(canonical body, min canonical slot, sorted canonical slots)` | rejected: body-primary is not *more* invariant than min-slot (already total), costs extra byte churn, and `sorted canonical slots` cannot ever be reached |
| widen the fold so the group-theoretic argument reaches `AGUKOD`-type drift | the relating transposition `(0 1)` is **not** in SPL's rotation group, so no amount of group enumeration reaches it; pinned by `test_canonical_slots.py::TestResidualClassIsOutOfReachByDesign::test_the_relating_transposition_is_not_a_proper_rotation_of_the_square` |
| fold same-symmetry-class, same-colour donors **within one fragment** (the known fix for the 23 `automorphism` cases) | **specified but deliberately NOT implemented.** It widens the fold beyond the geometry's own symmetry — the one boundary the whole design rests on — so it needs its own guards (a Δ/Λ tris-chelate divergence fixture at minimum) and a product call, not a unilateral extension |
| `CanonicalRankAtoms(breakTies=True)` as the instrument for "are these two donors interchangeable?" | wrong instrument: it settles symmetry ties on the *input index*, the very dependence under measurement. `tools/slot_drift_mechanism.py` uses `breakTies=False` classes instead |
| assuming `slot_renumber` rising meant a regression | per-molecule accounting: 8 + 11 molecules moved `rdkit_canonical → slot_renumber`, **0** the other way, **0** stable→drift. Reclassification |
| an angle-grid Euler search for the rotation group | `[0,90,120,180,240,270]` cannot express PBP's 72° five-fold; found 2 of 10. Pinned against return by `TestAlignerSymmetriesAreUnified::test_pbp_five_fold_is_now_reachable`, which asserts the C5 generator `(0,1,3,4,5,6,2)` is present |
| running the suite with `OIN_CANONICAL_SLOTS` exported before promotion | fails exactly `test_regression_stability.py::{test_cisplatin, test_transplatin, test_fac_irppy3}` and nothing else — stale hardcoded goldens for a lever that deliberately relabels, **not** an isomer merge. Left as-is by the lane on purpose; re-pinned at promotion |
| tests spelling "lever off" by **deleting** the env var | after promotion that means ON. 17 failures in v0.4.5 (6 more in v0.4.6). `test_canonical_slot_invariance.py` now writes `OIN_CANONICAL_SLOTS="0"` explicitly, with the reason in a comment |

## Where it landed

- **Lever:** `OIN_CANONICAL_SLOTS`, **default-ON** since v0.4.5 (`oin/levers.py::_DEFAULT_ON`).
  `OIN_STABLE_METAL_AC` (this branch's AC valence-capping order fix, commit `8bf9df61`) was
  promoted alongside it; its veto instrument `tools/geometry_tag_shift.py --n 300` reported
  **298 molecules: 0 string changes, 0 `[M_XXX]` changes, 0 coordination-number changes**.
- **Code:** `src/oinsmiles/oin/canonical_slots.py`; `src/oinsmiles/oin/compare.py`
  (`_polyhedron_signature` now a wrapper, imports from `canonical_slots`);
  `src/oinsmiles/utils/xyz2mol.py` (`_canonical_pivot` line 1038, `_canonical_z_sign` line 1070,
  their call sites in `_align_to_pai` at 1195 / 1232, `get_input_order_key` at 1834, post-pass at
  2033); `src/oinsmiles/utils/oin_aligner.py::_brute_force_symmetries` (delegates to
  `derive_rotation_group`).
- **Guards:** `tests/unit/test_canonical_slots.py` — `TestVertexTableIsSingleSourceOfTruth`,
  `TestRotationGroupIsAGroup` (`test_expected_orders`, `test_closed_under_composition_and_inverses`,
  `test_spanning_geometries_admit_only_proper_rotations`, `test_octahedron_is_24_not_48`,
  `test_unknown_geometry_is_none`), `TestAlignerSymmetriesAreUnified`
  (`test_aligner_agrees_with_the_derived_group_everywhere`,
  `test_pbp_five_fold_is_now_reachable`), `TestLexMinSignature`
  (`test_fac_and_mer_do_not_collapse`, `test_rotated_relabeling_of_one_isomer_collapses`,
  `test_tie_break_is_deterministic_and_lex_smallest`), `TestCanonicalSlotPermutation`,
  `TestCanonicalizeOinSlots` (`test_every_rotation_of_one_labeling_gives_one_string`,
  `test_winding_character_is_preserved_verbatim`, `test_axial_suffix_is_carried_through`,
  `test_metal_fragment_stays_first`, `test_emitted_labeling_achieves_the_keys_lex_min`,
  `test_comparison_key_is_unchanged_by_the_post_pass`), `TestCanonicalSlotMap`
  (`test_per_donor_ambiguity_is_exactly_the_genuine_automorphism`),
  `TestResidualClassIsOutOfReachByDesign` (4 tests, to be **inverted not deleted** when a fix
  lands). End-to-end: `tests/unit/test_canonical_slot_invariance.py` —
  `TestCanonicalSlotsAreInvariant::{test_byte_identical_under_renumbering_and_rotation,
  test_renumbering_alone_is_the_hard_case}` over `CisPlatin.xyz`, `Cis-PtCl2(en).xyz`, `FeCO5.xyz`,
  `PtMeNH3ClBr-Cis.xyz`; `TestLeverOffIsByteIdentical::test_goldens`;
  `TestLeverOnGoldens::{test_moved_goldens_are_relabelings_not_new_isomers,
  test_already_canonical_fixtures_do_not_churn}`; `TestNoOverFolding` (3 tests). Plus the re-pins in
  `test_regression_stability.py` and
  `test_stereo_roundtrip_diagnostics.py::TestEtaRingCanonicalization::test_non_eta_fragment_order_is_inert`.
- **Commits.** Branch `swimlane/v045-lane2`, tip `840eab84`, merged into `main`: `cc0758d7`
  (extract `oin/canonical_slots`, return the lex-min permutation), `51770c93` (merge Lane 1),
  `3722b18e` (WIP post-pass — the agent halted mid-task, `"DUDREA still drifts"`), `8bf9df61`
  (`OIN_STABLE_METAL_AC`), `12898c5a` (measurement + "the plan's premise is insufficient"),
  `840eab84` (pin exactly which goldens move, and that each is key-identical). Promotion
  `1450b5ce`, release `0d165845`.
- **Delivered, measured:** +12 byte-stable molecules (88 → 100 of ~150), 0 regressions, 0
  isomer-level change, `fragment_reorder` and `winding_star_drift` held at 0, TD-005 closed.

## Open questions / for the next agent

1. **The 23 `automorphism` residuals.** The fix is specified: fold, in addition to the rotation
   group, permutations exchanging donors **within one fragment** that are (a) in the same
   `breakTies=False` symmetry class and (b) the same colour. Scoped that narrowly it cannot touch
   metal Δ/Λ (which lives in the arrangement *across* fragments) nor the 7 `DISTINCT_donors` cases
   (different classes → never folded). **Next measurement before implementing:** build a Δ/Λ
   tris-chelate divergence fixture (`tests/fixtures/ZUMNEC.xyz` exists from Lane 7) and confirm the
   two enantiomers stay distinct *raw and at key level* with the widened fold; then invert
   `TestResidualClassIsOutOfReachByDesign`.
2. **The 7 `DISTINCT_donors` molecules are a soundness bug, not a canonicality one** —
   `IMACOO_comp_0`, `RUBTIS_comp_0`, `VEXHIR_comp_0`, `ZACFER_comp_0`, `HAVGIW_comp_0`,
   `KISQAG_comp_0`, `ZOSNUS_comp_0`. One of each pair's two strings is *wrong*. Now owned by Lane 9
   (wrong-donor). **Next measurement:** for each, determine which of the two labelings matches the
   3D donor→vertex assignment, and whether the error is in the fit or in the marker placement.
3. **Metal `@` tags.** 0/150 molecules emit one today, so nothing is stale — but Lane 5's Δ/Λ
   descriptor (`OIN_EMIT_METAL_CONFIG`, held off) **must** be derived from the canonical arrangement
   via `canonical_slot_map`, never from pre-relabel slots. If a metal `@SPn`/`@OHn` tag is ever
   emitted, re-run the 150-molecule census first.
4. **The PBP unification is unconditional.** If a future release wants strict corpus-wide
   byte-identity with levers OFF, `_brute_force_symmetries`' delegation is the one line to gate;
   its measured cost is `IMACOO_comp_0` alone (1/149).
5. **Out of scope and confirmed to survive this lever:** bond-order/aromaticity perception
   order-dependence (`CEBVIR_comp_0`, Lane 1) and the chiral-tag parity class (`FEQFIS_comp_0`,
   Lane 8) — see `docs/RENUMBERING_INSTABILITY_v0.4.5.md`. `DUDREA_comp_0` is closed by
   `OIN_STABLE_METAL_AC`; **do not reopen it as a slot problem.**
