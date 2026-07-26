# The `DISTINCT_donors` class — measured, and it is not a soundness class

> **v0.4.5 Lane 9.** Lane 2 found 32 residual `slot_renumber` pairs its canonical-slot
> post-pass could not close, classified them at atom level, and split them **23 benign
> automorphism / 7 `DISTINCT_donors`** — the 7 described as *"the slots land on inequivalent
> donor atoms, so one of the two strings is WRONG. A soundness defect, not a canonicality
> one."* Lane 9 existed to find out which of the two strings was right.
>
> **Answer: for all 7, both are right.** There is no wrong-donor defect. The `DISTINCT_donors`
> verdict is a false positive of the classifier, and the 7 belong to the same determinism class
> as the 23. **No encoder code changed in this lane** — a fix that froze one of two equally
> faithful labelings would have been a fix to a non-defect.

---

## 1. Why the geometry cannot prefer one labeling over the other

A donor's slot integer comes from a Kabsch fit of the donor direction vectors onto an
idealized template (`OINDiscreteAligner._map_to_template`). That fit is **exactly degenerate
over the polyhedron's proper-rotation group**, by construction and not by accident: if `g`
permutes the template vertices as a rotation `R_g`, then

```
template[p ∘ g] == R_g · template[p]
```

so aligning the same donors to `p` and to `p∘g` are the same optimization up to a rotation of
the answer, and the optimal residual is *identical*. The fit therefore determines the
donor→vertex map only **up to that group**. Which member of the tied coset is emitted is
settled by `_map_to_template`'s keep-test:

```python
if rmsd < best_rmsd:      # strict <  ->  first-seen wins among exact ties
```

`candidates` arrives in lexicographic order over `itertools.permutations` of the donor
enumeration order, and that order is **atom order**. So renumbering the input file re-orders
the donors, re-orders the permutation enumeration, and a different — equally optimal — member
of the same coset wins.

**Consequence.** Two labelings related by a proper rotation are equally faithful to the 3D
structure. A proper rotation is a change of reference frame: it preserves every isomer-level
relation (cis/trans, fac/mer — which vertices are adjacent vs opposite) and it preserves every
reflection-odd descriptor (metal Δ/Λ, eta winding, an axial sign), because it is not a
reflection. Only a relabeling that is **not** in the group could be a soundness defect.

That is the test Lane 9 applied.

## 2. Ground truth, per molecule

Instrument: `tools/wrong_donor_groundtruth.py`. It intercepts `_align_to_pai` (for
`__origIdx`, the only surviving link from a mol index back to an xyz line — `get_tmc_mol`
rebuilds the molecule metal-fragment-first, so a mol index is *not* an xyz index),
`_reduce_hapticity` (to attach each haptic-reduced donor's global indices), and
`_map_to_template` (for the selected fit). It then re-scores **every** permutation with the
same `scipy Rotation.align_vectors` call the encoder makes, so the degeneracy figures cannot
be an artifact of `_candidate_permutations`' `1e-6` prefilter.

Levers = Lane 2's **stacked** arm, i.e. the arm the classification was made in:
`OIN_CANONICAL_BODY=1 OIN_CANONICAL_PERCEPTION=1 OIN_CANONICAL_SLOTS=1 OIN_STABLE_METAL_AC=1`.
Slots are reported as **geometric fit vertices** (`GT_RAW=1`) — that is what "which donor sits
at which coordination-template vertex" means. The encoder applies two *further* relabelings
afterwards (`_permute_and_serialize`'s lex-max and the Lane 2 post-pass); both are group
elements, so conjugating by them leaves group membership unchanged and the verdict intact.

| molecule | geo | best rssd | perms tied to 1e-9 | \|rot group\| | nearest NON-tied (gap) | π = base vertex → variant vertex, per physical donor | which string is right? |
|---|---|---:|---:|---:|---:|---|---|
| `IMACOO_comp_0` | PBP | 5.108894e-01 | **2** | 10 | 5.108895e-01 (**6.8e-08**) | **identity** | both — fit did not move |
| `RUBTIS_comp_0` | SPL | 1.753664e-01 | 8 | 8 | 1.870683e+00 (1.70) | **identity** | both — fit did not move |
| `VEXHIR_comp_0` | SPL | 6.747391e-02 | 8 | 8 | 1.966972e+00 (1.90) | **identity** | both — fit did not move |
| `ZACFER_comp_0` | TET | 4.147324e-01 | 12 | 12 | 2.205179e+00 (1.79) | **identity** | both — fit did not move |
| `HAVGIW_comp_0` | OCT | 1.075096e-01 | 24 | 24 | 1.945778e+00 (1.84) | `(2,3,4,5,0,1)` ∈ O | both — proper rotation |
| `KISQAG_comp_0` | OCT | 1.913860e-01 | 24 | 24 | 1.911687e+00 (1.72) | `(4,5,0,1,2,3)` ∈ O | both — proper rotation |
| `ZOSNUS_comp_0` | OCT | 2.209427e-01 | 24 | 24 | 1.886258e+00 (1.67) | `(5,4,2,3,0,1)`, `(2,3,5,4,1,0)` ∈ O | both — proper rotation |

`|Δ rssd|` between the two competing labelings, over all 8 base/variant pairs:
`0.0 · 1.1e-15 · 2.5e-15 · 3.4e-15 · 4.0e-15 · 6.6e-15 · 8.3e-15 · 1.2e-14`. **Exactly tied.**

For 6 of 7 the exactly-tied set is *precisely* the rotation group's orbit of the winner
(8/8, 8/8, 12/12, 24/24, 24/24, 24/24) and the nearest **non**-tied permutation is
**1.7–1.9 away**. So the coset is isolated by three orders of magnitude and the choice inside
it is pure tie-break. Verdict tally: `identity 4, in_rotation_group 4`, **`NOT_A_ROTATION` 0**.

**These are the same pairs Lane 2 classified, not lookalikes.** The renumberings here are the
instrument's own (seed 1234), so it is worth checking that they land on the same strings. They
do: for all 7 the base string is byte-identical to the probe's, **every** variant string judged
is byte-identical to one the probe recorded, and the coverage is complete — `|mine| == |probe|`
for all 7 (1 variant each, 2 for `ZOSNUS`). Verified against
`results-v0.4.5-lane2/stacked/canonicality_probe_all.json`.

**The donor→atom translation is self-checking.** Donors are matched across encodings through
`__origIdx` and then back through the renumbering permutation. If that translation were wrong,
the two donor *sets* would not coincide — they would be arbitrary atoms, including hydrogens,
which is exactly what an earlier buggy version produced. Measured: `set(base) == set(variant)`
for all 7. And in all 7 the donor count equals the vertex count (PBP 7/7, SPL 4/4, TET 4/4,
OCT 6/6), so π constrains **every** vertex and the group-membership test is exact rather than
permissive.

**Independent corroboration.** `canonical_roundtrip_key` already agrees: `key_equal` is `True`
for **14/14** drifted variants across the 7, and none of the 7 appears in the stacked arm's
`key_broken` list. Had any pair been a different isomer, the comparison key would have said so.

**The torsion oracle was not needed and was deliberately not used.** Every observed π is a
*proper* rotation, so no reflection was ever folded and no chirality claim is at stake. The
Y1/Lane 7 hazard (a rigid oracle calling an achiral molecule chiral) cannot arise from a test
whose whole content is "is this permutation in the proper-rotation subgroup".

## 3. Why the wrong assignment "wins" — two sub-mechanisms, both order-dependent

The lane's premise was that `_map_to_template` picks the wrong permutation. It is more
interesting than that: **for 4 of the 7 the geometric fit does not move at all.**

| sub-mechanism | molecules | where the drift enters |
|---|---|---|
| the fit picks a different tied group element | `HAVGIW`, `KISQAG`, `ZOSNUS` | `_map_to_template`'s strict `<` over 24 exactly-tied permutations, lex-first = atom order |
| the fit picks the **same** element (π = identity) | `IMACOO`, `RUBTIS`, `VEXHIR`, `ZACFER` | **downstream** of the fit |

For the second group the donor→vertex map is bit-identical and the drift is introduced by
`_permute_and_serialize`'s lex-max over the same rotation group, its "Homogeneous Sorting"
assignment of same-`chem_id` item sets to ranks, and which of two equivalent donors the
fragment body renders first. Captured pre-post-pass strings show it directly:

```
RUBTIS base pre : [Rh_SPL].CC(=N{0}…)…n{1}1.[CH]{2}1=[CH]{2>}CC[CH]{3>}=[CH]{3}CC1
RUBTIS got  pre : [Rh_SPL].CC(=N{0}…)…n{1}1.[CH]{3}1=[CH]{3>}CC[CH]{2>}=[CH]{2}CC1
                                              ^^ the chelate is IDENTICAL; only which COD
                                                 arm is written first changed
```

Fit vertices for both: vertex 2 = atoms `(56, 58)`, vertex 3 = atoms `(46, 48)`. So the
physical assignment is unchanged and what moved is **which of the two symmetry-equivalent COD
alkene arms RDKit writes first** — `_smilesAtomOutputOrder`, i.e.
`CanonicalRankAtoms(breakTies=True)` settling a tie between equivalent atoms on the input
index. **That is exactly the root cause Lane 2 named for the 23.** RUBTIS is in the 23's class.

Same shape for the others in that group: `VEXHIR`'s two acac oxygens (the `=O` position holds
atom 138 in one encoding and atom 134 in the other), `IMACOO`'s two β-diketiminate nitrogens
(amide vs imine), `ZACFER`'s two Cp\* rings.

**A note on three of the seven that matters for anyone re-reading the classifier's output.**
In `VEXHIR`/`KISQAG` (acetylacetonate), `IMACOO` (β-diketiminate) and `HAVGIW` (dipyrrin), the
two donors the classifier called "inequivalent" are inequivalent **only in the localized
Kekulé form the encoder happens to have written** — one O is a ketone and one an enol, one N
an amide and one an imine. The ligands are constitutionally symmetric and the donors are a
resonance pair. `CanonicalRankAtoms` is being asked about a localization, not about the
molecule.

## 4. Why the classifier said `DISTINCT_donors`

`tools/slot_drift_mechanism.py::atom_verdict` asks a **per-fragment, positionally-zipped**
question about a **global** relabeling. Three independent failure modes, each observed:

1. **A rotation acts on every fragment at once.** `RUBTIS`'s chelate donors really are
   inequivalent (imine N vs pyridyl N), so locally the swap looks illegitimate — but the full
   relabeling is `(0 1)(2 3)`, a C₂ of the square, because the COD's two arms swapped as well.
   A per-fragment test can never see the composition.
2. **`zip(frags_b, frags_g)` assumes fragment order is stable.** It is not: the post-pass
   re-derives fragment order *from the slot integers*, so any molecule whose slots moved can
   have its fragments reordered — and a complex with two copies of one ligand then gets ligand
   A compared against ligand B.
3. **Haptic donors were keyed on the slot marker's *first occurrence*.** An η⁵ ring stamps its
   slot on every ring atom; comparing whichever ring position was written first makes two
   genuinely interchangeable Cp\* rings present a methyl-bearing carbon against a
   silyl-bearing one, so they land in different symmetry classes for no chemical reason.
   `ZACFER` is exactly this.

A fourth, structural limitation surfaced while measuring and is **not** fixed here: the
verdict logic asks "does slot *s* sit on an equivalent atom in both strings?", which presumes
slot *s* stays in the same fragment. In `ZOSNUS` slots 4 and 5 move *between* the two dppe
ligands, so no per-fragment question can resolve it. Closing that needs the group-membership
test, i.e. the 3D instrument or a string-level implementation of it.

## 5. What changed in this lane

**No encoder change. No new lever.** Levers-OFF bytes are unchanged because no encoder file
was touched; the fix a soundness class would have needed does not apply to a class that is
sound. Changed instead:

| file | change |
|---|---|
| `tools/wrong_donor_groundtruth.py` | new. The geometric ground-truth instrument, with the degeneracy argument in its docstring so the next reader does not re-derive it. |
| `tools/slot_drift_mechanism.py` | `DISTINCT_donors` → `distinct_donors_LOCAL`, with a warning that it is a local string-only heuristic and **cannot** establish that either string is wrong; fragments paired by **body text** instead of position; haptic donors compared on the **whole** atom set instead of the first occurrence. |
| `tests/unit/test_wrong_donor_verdict.py` | new. Pins the two mechanisms fixed, using the real strings from the actual 7, plus the group-theoretic fact that makes the per-fragment test misleading. |

Measured effect on the stacked arm, per molecule, old classifier vs new:

```
pairs classified: 32   unchanged verdict: 25
OLD automorphism: 23   NEW automorphism: 24
moved OUT of automorphism (regression if non-empty): []          <- the 23 did not move
moved INTO automorphism: ZACFER_comp_0                            <- haptic first-occurrence fix
remaining flagged (rename only): IMACOO RUBTIS VEXHIR HAVGIW KISQAG ZOSNUS
```

The fragment-body pairing fix changed **0** verdicts on this corpus — it is defensive, not
load-bearing, and is reported as such rather than claimed as a win.

## 6. What to do instead — and the blocker that is now removed

The real defect behind all 30 (23 + 7) is one thing: **the choice among group-related
labelings is settled on atom order**, in three places (`_map_to_template`'s strict `<`,
`_permute_and_serialize`'s lex-max, and the fragment body's atom output order).

Lane 2's `docs/CANONICAL_SLOTS_v0.4.5.md` §7a specifies the fold that closes this and declines
it for three stated reasons. Lane 9's measurement removes one of them outright:

- ~~*"it would over-fold your 7"*~~ — **refuted.** All 7 are already related by a proper
  rotation, which the post-pass folds over anyway. There is nothing there to over-fold.
- *"it folds past the geometry's own symmetry"* — still true and still the real risk, but note
  the geometry's own symmetry is `|G|`-fold **exactly**, and every π measured here is inside it.
- *"needs a Δ/Λ guard and a product call"* — unchanged. Still not an agent's call to make.

A fix aimed only at the Kabsch tie-break would close at most **3** of the 7; the other 4 need
the downstream tie-breaks made invariant too. Do **not** ship a lever that merely makes the 7
deterministic on the premise that one of them was wrong.

## 7. Loose ends, honestly

- **PBP template rounding.** `IMACOO` is the only case whose exactly-tied set (2) is not the
  full group orbit (10), and its nearest non-tied permutation is **6.8e-08** away.
  `GEOMETRY_VERTICES["PBP"]` is stored to 7 decimals, so the derived C₅ is only approximately
  a symmetry of the *stored* template and the orbit splits at ~1e-7. For PBP the tie is decided
  by template rounding noise rather than by the molecule. Separate order-sensitivity, unchased.
  (`IMACOO` is also the single molecule Lane 2's ungated rotation-group unification moved.)
- **2 of the 32 residuals are `unparsable`** (borane-cluster class) and have never been
  classified at atom level by anyone. Untouched here.
- `postpass_BUG_diverges` is `4` on the levers-OFF arm and `3` on the Lane-1-only arm, against
  Lane 2's *"Measured 0"*. Those arms' strings were not produced by the post-pass, so this is
  the post-pass failing to fold strings it never emitted, not a new defect — but the "0" claim
  holds only for the `on`/`stacked` arms and should be read that way.
