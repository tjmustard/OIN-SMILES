# Lane 9 — the seven molecules with slots on inequivalent donors

**What this lane was for:** a user-requested **correctness** investigation into the 7 residual
`slot_renumber` pairs Lane 2 had classified as `DISTINCT_donors` — *"the slots land on inequivalent
donor atoms, so one of the two strings is WRONG. A soundness defect, not a canonicality one."*
Lane 9 existed to find out **which** of the two strings was right.

**Verdict: for all 7, both are.** `DISTINCT_donors` is **not a soundness class**. It is a false
positive of a string-only classifier, and the 7 belong to the **same determinism class as the 23**
benign pairs. **No encoder code changed in this lane** — a fix would have frozen one of two equally
faithful labelings, i.e. fixed a non-defect.

The seven: **`IMACOO`, `RUBTIS`, `VEXHIR`, `ZACFER`, `HAVGIW`, `KISQAG`, `ZOSNUS`** (each
`_comp_0`).

Primary sources: `docs/WRONG_DONOR_v0.4.5.md`, `tests/unit/test_wrong_donor_verdict.py`,
`tools/wrong_donor_groundtruth.py`, `docs/V045_STATUS_2026-07-25.md` §"RESOLVED — the '7
wrong-donor molecules' were NOT a defect".

**This is a clean negative result, and it is reported as one.** The suspicion was reasonable, the
measurement dismissed it, and the value of the lane is (a) that no code was written to fix a
non-defect, (b) that the classifier which raised the alarm no longer makes a claim it cannot
support, and (c) that one of the three reasons Lane 2 gave for declining a wider fold is now
refuted, which changes the shape of the *real* fix.

---

## ELI5

The 1D string labels each atom bound to the metal with a slot number — "this nitrogen is at
position 0, this chlorine at position 3" — and the numbers come from fitting the real 3D shape onto
an idealised template (an octahedron, a square, and so on). Lane 2 noticed that if you shuffle the
lines of the input file, some molecules come back with the slot numbers landing on *different*
atoms, and flagged 7 cases where the two atoms that swapped were chemically different from each
other — which would mean the string was lying about which atom is bound where. That would be
serious. Lane 9 checked it against the coordinates and found the alarm was false: the two labelings
always describe the **same** arrangement seen from a **rotated viewpoint**, and rotating your
viewpoint does not change a molecule. Numbering an octahedron's corners differently after turning
it a quarter-turn is not a mistake; it is the same octahedron. The classifier could not see this
because it examined one ligand at a time, and a rotation acts on **all** the ligands at once.

---

## The work, visually

```
PART 1 — WHY THE GEOMETRY CANNOT PREFER ONE LABELING OVER THE OTHER

A donor's slot integer comes from a Kabsch fit of the donor DIRECTION VECTORS onto an
idealized template (utils/oin_aligner.py :: OINDiscreteAligner._map_to_template).

That fit is EXACTLY DEGENERATE over the polyhedron's PROPER-ROTATION GROUP,
by construction and not by accident:

        if g permutes the template vertices as a rotation R_g, then

                     template[p ∘ g]  ==  R_g · template[p]

        so aligning the same donors to  p  and to  p∘g  is THE SAME OPTIMIZATION
        up to a rotation of the answer, and the optimal residual is IDENTICAL.

  ⇒ the fit determines donor → vertex ONLY UP TO THAT GROUP.

  which tied member is emitted is settled here:

        if rmsd < best_rmsd:      # strict  <   ⇒ FIRST-SEEN wins among exact ties
                                  #             ⇒ lex over itertools.permutations
                                  #             ⇒ lex over the donor enumeration order
                                  #             ⇒ ATOM ORDER

  renumber the input file  →  donors re-order  →  permutation enumeration re-orders
                           →  a DIFFERENT, EQUALLY OPTIMAL member of the same coset wins.


   base labeling                        variant labeling
        1                                    3
        │                                    │
   0 ──(M)── 2      ==  rotate 90° ==    2 ──(M)── 0
        │                                    │
        3                                    1

   A PROPER ROTATION IS A CHANGE OF REFERENCE FRAME.  It preserves:
        · every isomer-level relation   (cis/trans, fac/mer — which vertices are
          adjacent vs opposite)
        · every reflection-ODD descriptor (metal Δ/Λ, eta winding, an axial sign),
          BECAUSE IT IS NOT A REFLECTION
   ⇒ ONLY a relabeling NOT in the group could be a soundness defect.
     That is the test Lane 9 applied.


PART 2 — WHY THE CLASSIFIER SAID `DISTINCT_donors`
          (tools/slot_drift_mechanism.py :: atom_verdict)

It asks a PER-FRAGMENT, POSITIONALLY-ZIPPED question about a GLOBAL relabeling.

   RUBTIS_comp_0, Rh square planar, COD + pyridyl-imine chelate:

     what the classifier SEES (one fragment at a time)      what ACTUALLY happened
     ──────────────────────────────────────────────────     ───────────────────────────
     chelate:  N{2} ↔ n{3}                                  chelate:  (2 3)
        imine N vs pyridyl n — GENUINELY inequivalent        COD arms: (0 1)   ← equivalent,
        ⇒ "(2 3)" is NOT a rotation of the square                        invisible per-fragment
        ⇒ FLAG: DISTINCT_donors  ✗                          composite: (0 1)(2 3)  ∈ D4  ✓
                                                             ⇒ a PROPER ROTATION

   three independent mechanical failure modes, each observed:
     (1) a rotation acts on EVERY fragment at once      → composition is invisible per-fragment
     (2) zip(frags_b, frags_g) assumes fragment ORDER is stable — the post-pass RE-DERIVES
         fragment order FROM the slot integers, so two copies of one ligand get compared
         against each other
     (3) haptic donors were keyed on the slot marker's FIRST OCCURRENCE — an η5 ring stamps
         its slot on every ring atom, so two interchangeable Cp* rings presented a
         methyl-bearing carbon against a silyl-bearing one   (ZACFER is exactly this)
   + (4) a STRUCTURAL limit, NOT fixed: the verdict presumes slot s stays in the same
         fragment.  In ZOSNUS slots 4 and 5 move BETWEEN the two dppe ligands, so NO
         per-fragment question can resolve it.


PART 3 — WHAT WAS MEASURED  (tools/wrong_donor_groundtruth.py)

  xyz ──► encode base ──┐
      └─► renumber ─────┴─► encode variant
                                │
        intercept 3 encoder call sites (each patch WRAPS the original, returns its result):
              _align_to_pai      → __origIdx, the ONLY surviving link from a mol index
                                   back to an xyz line (get_tmc_mol rebuilds metal-fragment
                                   -first, so a mol index is NOT an xyz index)
              _reduce_hapticity  → each haptic-reduced donor's global indices
              _map_to_template   → the fit the encoder actually SELECTED
                                │
        re-score EVERY permutation with the SAME scipy Rotation.align_vectors the encoder
        uses  (bypasses _candidate_permutations, so the degeneracy figures cannot be an
        artifact of its 1e-6 prefilter)
                                │
        π  =  (base vertex → variant vertex) per PHYSICAL donor
                                │
        ┌───────────────────────┴─────────────────────────────────────────┐
        │ π ∈ rotation group  → determinism defect. BOTH strings faithful │
        │ π == identity       → the fit did not move; drift is DOWNSTREAM │
        │ π ∉ rotation group  → a GENUINE soundness defect: investigate   │
        └─────────────────────────────────────────────────────────────────┘

                 RESULT:  identity 4 · in_rotation_group 4 · NOT_A_ROTATION 0

LEGEND
  ✓ / ✗          the correct reading / the classifier's false reading
  π              the permutation taking each physical donor's base vertex to its variant vertex
  rssd           the Kabsch fit's residual (scipy's `align_vectors` root-sum-square deviation)
  D4             the proper-rotation group of the square (SPL); |D4| = 8
  coset          the set of labelings tied at the optimal residual
  the 23 / the 7 Lane 2's split of its 32 residual `slot_renumber` pairs
  §7a            docs/CANONICAL_SLOTS_v0.4.5.md §7a — the wider fold Lane 2 specified and declined
```

---

## Initial assumptions and hypothesis

**The suspicion was reasonable, and it is worth being explicit about why**, because "clean
negative" can otherwise read as "someone wasted a lane."

1. **Lane 2 had already done the honest thing.** Its canonical-slot post-pass closed most slot
   drift; it then classified the 32 residual pairs *at atom level* rather than declaring victory,
   and it **separated** 23 pairs it could explain (benign automorphism) from 7 it could not. That
   is the behaviour you want.
2. **The 7's description was a real, checkable claim.** "The slots land on inequivalent donor
   atoms" is not vague — for `RUBTIS` the chelate's two donors genuinely *are* inequivalent (imine
   N vs pyridyl n). If slot 2 sits on the imine N in one string and on the pyridyl n in the other,
   then on the face of it one string is wrong about which atom is bound where. That is a soundness
   defect and would have been the most serious finding of the release.
3. **The hypothesis had a concrete suspect.** `_map_to_template` picks the permutation, and its
   keep-test is a strict `<`, which among exact ties keeps whichever candidate was enumerated
   first. Atom order decides. That is a textbook non-determinism, and "the fit picks the wrong
   permutation" was the natural reading.
4. **The prior was also that a fix existed.** Lane 2's `docs/CANONICAL_SLOTS_v0.4.5.md` §7a already
   specified a wider fold and declined it for three stated reasons, one of which was **"it would
   over-fold your 7"**. So the 7 were simultaneously the reason not to take the fix and the reason
   a fix was needed.

The hypothesis Lane 9 set out to test was therefore: **`_map_to_template` selects a labeling that
places a slot on the wrong physical donor, and for each of the 7 exactly one of the two strings is
faithful to the 3D structure.**

The design decision that made it answerable: **settle it from the coordinates, not from the
strings.** `slot_drift_mechanism.py` has only the strings, and "which physical atom is bound where"
is not a question strings can answer.

---

## What was actually found

### Method and configuration, stated so the numbers can be reproduced

Instrument: `tools/wrong_donor_groundtruth.py`. Levers set to **Lane 2's *stacked* arm** — the arm
the classification was made in:

```
OIN_CANONICAL_BODY=1 OIN_CANONICAL_PERCEPTION=1 OIN_CANONICAL_SLOTS=1 OIN_STABLE_METAL_AC=1
```

Slots reported as **geometric fit vertices** (`GT_RAW=1`) — that is what *"which donor sits at
which coordination-template vertex"* means. The encoder applies two **further** relabelings
afterwards (`_permute_and_serialize`'s lex-max, and the Lane 2 post-pass); **both are group
elements**, so conjugating by them leaves group membership unchanged and the verdict intact.
Renumberings are the instrument's own, **seed 1234**.

### CONFIRMED — ground truth, per molecule: 0 soundness defects

| molecule | geo | best rssd | perms tied to 1e-9 | \|rot group\| | nearest NON-tied (gap) | π = base vertex → variant vertex, per physical donor | which string is right? |
|---|---|---:|---:|---:|---:|---|---|
| `IMACOO_comp_0` | PBP | 5.108894e-01 | **2** | 10 | 5.108895e-01 (**6.8e-08**) | **identity** | both — fit did not move |
| `RUBTIS_comp_0` | SPL | 1.753664e-01 | 8 | 8 | 1.870683e+00 (1.70) | **identity** | both — fit did not move |
| `VEXHIR_comp_0` | SPL | 6.747391e-02 | 8 | 8 | 1.966972e+00 (1.90) | **identity** | both — fit did not move |
| `ZACFER_comp_0` | TET | 4.147324e-01 | 12 | 12 | 2.205179e+00 (1.79) | **identity** | both — fit did not move |
| `HAVGIW_comp_0` | OCT | 1.075096e-01 | 24 | 24 | 1.945778e+00 (1.84) | `(2,3,4,5,0,1)` ∈ O | both — proper rotation |
| `KISQAG_comp_0` | OCT | 1.913860e-01 | 24 | 24 | 1.911687e+00 (1.72) | `(4,5,0,1,2,3)` ∈ O | both — proper rotation |
| `ZOSNUS_comp_0` | OCT | 2.209427e-01 | 24 | 24 | 1.886258e+00 (1.67) | `(5,4,2,3,0,1)`, `(2,3,5,4,1,0)` ∈ O | both — proper rotation |

**Tally: `identity` 4, `in_rotation_group` 4, `NOT_A_ROTATION` 0.**

> Note on counting, because two different totals appear in the source doc and they count different
> things. There are **7 molecules** but **8 base/variant pairs judged** — `ZOSNUS` contributes two
> variants, which is why the tally sums to 8 and why two π values are listed for it. Separately,
> `canonical_roundtrip_key` agreement is reported over **14/14 drifted variants**, a larger set
> drawn from the probe's own renumberings rather than the 8 judged pairs.

### CONFIRMED — the two competing labelings are EXACTLY tied

`|Δ rssd|` between the two competing labelings, over all **8** base/variant pairs:

```
0.0 · 1.1e-15 · 2.5e-15 · 3.4e-15 · 4.0e-15 · 6.6e-15 · 8.3e-15 · 1.2e-14
```

For **6 of 7** the exactly-tied set is *precisely* the rotation group's orbit of the winner
(**8/8, 8/8, 12/12, 24/24, 24/24, 24/24**) and the nearest **non**-tied permutation is
**1.7–1.9 away**. So the coset is isolated by roughly **three orders of magnitude** and the choice
inside it is **pure tie-break**.

### CONFIRMED — these are the same pairs Lane 2 classified, not lookalikes

Because the renumberings are the instrument's own, this had to be checked rather than assumed. It
holds: for **all 7** the base string is byte-identical to the probe's, **every** variant string
judged is byte-identical to one the probe recorded, and coverage is **complete** — `|mine| ==
|probe|` for all 7 (1 variant each, 2 for `ZOSNUS`). Verified against
`results-v0.4.5-lane2/stacked/canonicality_probe_all.json`.

### CONFIRMED — the donor → atom translation is self-checking, and the group test is exact

If the `__origIdx`-based translation were wrong, the two donor **sets** would not coincide — they
would be arbitrary atoms, including hydrogens, **which is exactly what an earlier buggy version
produced**. Measured: `set(base) == set(variant)` for **all 7**. And in all 7 the **donor count
equals the vertex count** (PBP 7/7, SPL 4/4, TET 4/4, OCT 6/6), so π constrains **every** vertex
and the group-membership test is **exact rather than permissive**.

### CONFIRMED — independent corroboration from the key

`canonical_roundtrip_key` already agrees: **`key_equal` is `True` for 14/14** drifted variants
across the 7, and **none of the 7 appears in the stacked arm's `key_broken` list.** Had any pair
been a different isomer, the comparison key would have said so.

### CONFIRMED — a sub-mechanism nobody predicted: for 4 of 7 the fit does not move AT ALL

The lane's premise was that `_map_to_template` picks the wrong permutation. It is more interesting
than that.

| sub-mechanism | molecules | where the drift enters |
|---|---|---|
| the fit picks a different tied group element | `HAVGIW`, `KISQAG`, `ZOSNUS` | `_map_to_template`'s strict `<` over **24 exactly-tied** permutations, lex-first = atom order |
| the fit picks the **same** element (π = identity) | `IMACOO`, `RUBTIS`, `VEXHIR`, `ZACFER` | **downstream** of the fit |

For the second group the donor → vertex map is **bit-identical** and the drift is introduced by
`_permute_and_serialize`'s lex-max over the same rotation group, its "Homogeneous Sorting"
assignment of same-`chem_id` item sets to ranks, and **which of two equivalent donors the fragment
body renders first.** Captured pre-post-pass strings show it directly:

```
RUBTIS base pre : [Rh_SPL].CC(=N{0}…)…n{1}1.[CH]{2}1=[CH]{2>}CC[CH]{3>}=[CH]{3}CC1
RUBTIS got  pre : [Rh_SPL].CC(=N{0}…)…n{1}1.[CH]{3}1=[CH]{3>}CC[CH]{2>}=[CH]{2}CC1
                                              ^^ the chelate is IDENTICAL; only which COD
                                                 arm is written first changed
```

Fit vertices for both: **vertex 2 = atoms (56, 58), vertex 3 = atoms (46, 48).** So the physical
assignment is unchanged and what moved is **which of the two symmetry-equivalent COD alkene arms
RDKit writes first** — `_smilesAtomOutputOrder`, i.e. `CanonicalRankAtoms(breakTies=True)` settling
a tie between equivalent atoms **on the input index**. **That is exactly the root cause Lane 2
named for the 23.** `RUBTIS` is in the 23's class.

Same shape for the rest of that group: `VEXHIR`'s two acac oxygens (the `=O` position holds atom
138 in one encoding and atom 134 in the other), `IMACOO`'s two β-diketiminate nitrogens (amide vs
imine), `ZACFER`'s two Cp\* rings.

### CONFIRMED — three of the seven were never inequivalent in the first place

In `VEXHIR`/`KISQAG` (acetylacetonate), `IMACOO` (β-diketiminate) and `HAVGIW` (dipyrrin), the two
donors the classifier called "inequivalent" are inequivalent **only in the localized Kekulé form
the encoder happens to have written** — one O is a ketone and one an enol, one N an amide and one
an imine. **The ligands are constitutionally symmetric and the donors are a resonance pair.**
`CanonicalRankAtoms` is being asked about a *localization*, not about the molecule.

### CONFIRMED — the torsion oracle was NOT needed and was deliberately not used

Every observed π is a **proper** rotation, so no reflection was ever folded and **no chirality
claim is at stake**. The Y1/Lane 7 hazard — a rigid oracle calling an achiral molecule chiral —
cannot arise from a test whose whole content is *"is this permutation in the proper-rotation
subgroup?"*

### CONFIRMED — the classifier A/B, including the part that changed nothing

Measured effect on the stacked arm, per molecule, old classifier vs new:

```
pairs classified: 32   unchanged verdict: 25
OLD automorphism: 23   NEW automorphism: 24
moved OUT of automorphism (regression if non-empty): []          <- the 23 did not move
moved INTO automorphism: ZACFER_comp_0                            <- haptic first-occurrence fix
remaining flagged (rename only): IMACOO RUBTIS VEXHIR HAVGIW KISQAG ZOSNUS
```

The **fragment-body pairing fix changed 0 verdicts on this corpus** — reported as **defensive, not
load-bearing**, rather than claimed as a win.

---

## What was done

### Nothing, in `src/` — and that is the deliverable

**No encoder change. No new lever.** Levers-OFF bytes are unchanged **by construction** because no
file under `src/` was touched. The fix a soundness class would have needed does not apply to a class
that is sound, and shipping one would have frozen one of two equally faithful labelings.

### The instrument — `tools/wrong_donor_groundtruth.py` (new, 341 lines)

The degeneracy argument lives in the module docstring **so the next reader does not have to
re-derive it.** It patches three encoder call sites; **every patch wraps the original and returns
its result**, so nothing about encoder behaviour changes:

| patched | why |
|---|---|
| `_align_to_pai` | for `__origIdx` — the **only** surviving link from a mol index back to an xyz line, because `get_tmc_mol` rebuilds the molecule **metal-fragment-first**, so a mol index is **not** an xyz index |
| `_reduce_hapticity` | to attach each haptic-reduced donor's global indices |
| `_map_to_template` | for the fit the encoder actually **selected** |

It then **re-scores every permutation** with the same `scipy` `Rotation.align_vectors` call the
encoder makes, deliberately bypassing `_candidate_permutations` so the degeneracy figures cannot be
an artifact of its `1e-6` prefilter. `GT_RAW=1` reports the raw fit vertices; without it the Lane 2
post-pass map is composed in.

Usage:

```bash
PYTHONPATH=src OIN_CANONICAL_BODY=1 OIN_CANONICAL_PERCEPTION=1 \
OIN_CANONICAL_SLOTS=1 OIN_STABLE_METAL_AC=1 GT_RAW=1 \
    .venv/bin/python tools/wrong_donor_groundtruth.py [all|<molecule>] [n_renumberings]
```

**Rejected alternative:** answering the question from the strings, by extending
`slot_drift_mechanism.py`. Rejected on principle — *"which physical atom is bound where"* is a
claim about coordinates, and a tool that only has strings **cannot** settle it. That is the
lane's central methodological point, and the tool rename below is how it was enforced rather than
merely noted.

### The classifier repair — `tools/slot_drift_mechanism.py` (Lane 2's tool)

| change | why |
|---|---|
| `DISTINCT_donors` → **`distinct_donors_LOCAL`**, with an attached warning | it is a **local, string-only heuristic** and **cannot** establish that either string is wrong. The old name asserted a soundness defect the tool has no way to see |
| fragments paired by **body text** instead of by position (`_pair_fragments`) | the post-pass re-derives fragment order **from the slot integers**, so any molecule whose slots moved can have its fragments reordered — and a complex with two copies of one ligand then gets ligand A compared against ligand B |
| haptic donors compared on the **whole** atom set instead of the marker's first occurrence (`_slot_to_atoms`) | an η⁵ ring stamps its slot on **every** ring atom; comparing whichever position was written first made two genuinely interchangeable Cp\* rings present a methyl-bearing carbon against a silyl-bearing one. **`ZACFER` is exactly this** |

**Housekeeping decision recorded in the status doc:** Lane 9's edit to Lane 2's tool **stays on
Lane 9's branch**. Both merge into `main`, so moving commits between branches adds rebase risk for
no benefit, and **Lane 9 is the lane whose evidence justifies the change.**

### The guard — `tests/unit/test_wrong_donor_verdict.py` (new)

Uses the **real emitted strings** from the actual 7, so it needs **no dataset** (the dataset is
gitignored and absent from a worktree). Two full string pairs are inlined: `ZACFER_BASE`/
`ZACFER_GOT` (two equivalent Cp\* rings bridged by Si–Si on a TET Ti centre) and `RUBTIS_BASE`/
`RUBTIS_GOT` (COD + pyridyl-imine chelate on a SPL Rh centre).

| class / test | what it pins |
|---|---|
| `TestBothStringsAreTheSameMolecule::test_both_pairs_are_same_vcolor_identical` | the first-stage class must stay `same_vcolor_identical` — the premise of stage two |
| `TestHapticFirstOccurrenceFalsePositive::test_zacfer_is_an_automorphism` | `atom_verdict(ZACFER_BASE, ZACFER_GOT) == "automorphism"` |
| `::test_slot_markers_collect_every_occurrence` | a 5-atom eta group must contribute **5** atoms to its slot, not 1 |
| `TestFragmentsPairedByBodyNotPosition::test_reordered_fragments_still_pair_on_body` | same bodies in a different order must still pair correctly |
| `::test_absent_body_is_reported_not_mispaired` | `_pair_fragments(["N{0}CC"], ["O{0}CC"])` is `None` |
| `TestVerdictMakesNoSoundnessClaim::test_rubtis_verdict_is_the_local_heuristic_name` | the verdict is `distinct_donors_LOCAL` and **does not contain** `DISTINCT_donors` |
| `::test_module_does_not_reintroduce_the_soundness_claim` | greps `slot_drift_mechanism.py`'s own source for `return "DISTINCT_donors"` — **the name cannot come back** |
| `TestWhyAPerFragmentTestIsMisleading::test_chelate_swap_alone_is_not_a_rotation_but_the_composite_is` | `(0,1,3,2)` = `(2 3)` is **not** in SPL's group; `(1,0,3,2)` = `(0 1)(2 3)` **is**. The group-theoretic fact behind RUBTIS, pinned so the reasoning survives |
| `::test_spanning_geometries_admit_only_proper_rotations` | `det > 0` for every element of the OCT, TET and PBP groups, with a floor assertion of **≥ 46 permutations checked** (24 + 12 + 10) so the test cannot pass by silently checking nothing |

**Commit `c804044f` — `test(lane9): state the proper-rotation guard's scope instead of implying
it`** exists for that last test and is worth citing on its own. The test **skipped planar templates
silently**, which *reads* as an oversight. It is not: for a **rank-3** (spanning) vertex set the
realizing linear map is **unique**, so `det > 0` is a real assertion — but for a **planar** set
(SPL, TPL, LIN) the out-of-plane direction is **free**, every Gram-preserving permutation extends to
a proper 3D rotation, and the check would be **vacuous**. `(0 1)` on the square is the concrete
case: a C₂ about the in-plane axis through vertices 2 and 3, genuinely proper. The commit says so in
the docstring, hoists the basis search out of the per-permutation loop (it does not depend on the
permutation), and adds the ≥ 46 floor. It also states the test's own status honestly: it **mirrors**
`derive_rotation_group`'s `det > 0` filter, so it is a guard against that filter being removed, **not
an independent derivation** — and it is there *because the Y2 axial wave lost chirality exactly by
folding over something that turned out to be a reflection.*

---

## Dead ends and refutations

**The entire lane is a negative result.** It is presented as the finding, not as a shortfall: the
lane's product is a dismissed defect, a disarmed classifier, and a corrected input to the *real*
fix decision.

### REFUTED — `DISTINCT_donors` is a soundness class

`identity 4 / in_rotation_group 4 / NOT_A_ROTATION 0`, `|Δrssd| ≤ 1.2e-14`, `key_equal` True
14/14, donor count == vertex count for all 7. The 7 are in the same determinism class as the 23.
`docs/V045_STATUS_2026-07-25.md`'s section on the "7 wrong-donor molecules" is annotated
**"✅ RESOLVED — NOT a defect"** with an instruction to read `docs/WRONG_DONOR_v0.4.5.md` before
acting on anything under that heading, **because the wrong version was reported and acted on.**

### REFUTED — "the fix is in `_map_to_template`'s tie-break"

For **4 of the 7 the geometric fit does not move at all** (π = identity). **A fix aimed only at the
Kabsch tie-break would close at most 3 of the 7**; the other 4 need the downstream tie-breaks made
invariant too. And the corollary that matters most: **do not ship a lever that merely makes the 7
deterministic on the premise that one of them was wrong.** Determinism is a legitimate goal; the
premise is not.

### REFUTED — one of Lane 2's three reasons for declining the §7a wider fold

Lane 2's `docs/CANONICAL_SLOTS_v0.4.5.md` §7a specifies the fold that closes this seam and declines
it for three stated reasons. Lane 9 removes one **outright**:

- ~~*"it would over-fold your 7"*~~ — **REFUTED.** All 7 are **already** related by a proper
  rotation, which the post-pass folds over anyway. **There is nothing there to over-fold.**
- *"it folds past the geometry's own symmetry"* — **still true and still the real risk**, but note
  the geometry's own symmetry is `|G|`-fold **exactly**, and every π measured here is inside it.
- *"needs a Δ/Λ guard and a product call"* — **unchanged. Still not an agent's call to make.**

The recommendation in the status doc remains **not** to take the fold unilaterally.

### The real defect, named — and it is one sentence in three places

**The choice among group-related labelings is settled on atom order.** Three sites:

1. `_map_to_template`'s strict `<` (keeps the lex-first of exactly-tied candidates);
2. `_permute_and_serialize`'s lex-max over the same rotation group (plus its homogeneous-sorting
   rank assignment);
3. the fragment body's atom output order — `CanonicalRankAtoms(breakTies=True)` /
   `_smilesAtomOutputOrder` settling a tie between equivalent atoms on the input index.

That single sentence is the root cause behind **all 30** residuals (23 + 7), not just the 7.

### Loose ends, carried honestly

- **PBP template rounding.** `IMACOO` is the **only** case whose exactly-tied set (**2**) is not the
  full group orbit (**10**), and its nearest non-tied permutation is **6.8e-08** away.
  `GEOMETRY_VERTICES["PBP"]` is stored to **7 decimals**, so the derived C₅ is only *approximately*
  a symmetry of the **stored** template and the orbit splits at ~1e-7. **For PBP the tie is decided
  by template rounding noise rather than by the molecule.** A separate order-sensitivity, unchased.
  (`IMACOO` is also the single molecule Lane 2's ungated rotation-group unification moved.)
- **2 of the 32 residuals are `unparsable`** (borane-cluster class) and have **never** been
  classified at atom level by anyone. Untouched here.
- **`postpass_BUG_diverges` is 4 on the levers-OFF arm and 3 on the Lane-1-only arm**, against Lane
  2's *"Measured 0"*. Those arms' strings were **not produced by the post-pass**, so this is the
  post-pass failing to fold strings it never emitted, **not** a new defect — but **Lane 2's "0"
  claim holds only for the `on`/`stacked` arms and should be read that way.**
- **The fourth classifier failure mode is structural and NOT fixed.** `atom_verdict` asks *"does
  slot `s` sit on an equivalent atom in both strings?"*, which presumes slot `s` stays in the same
  fragment. **In `ZOSNUS` slots 4 and 5 move *between* the two dppe ligands**, so no per-fragment
  question can resolve it. Closing that needs the group-membership test — i.e. the 3D instrument, or
  a string-level implementation of it.

---

## Where it landed

**Branch `swimlane/v045-lane9`, tip `c804044f`. Fully merged** — `main` is 159 commits ahead and 0
behind, so `git log main..swimlane/v045-lane9` is empty. Merged via `14a761f6`
(→ `trial/v045-integration`) and `bbbfb3f8` (→ `release/v0.4.5`).

Commits, both of them:

| commit | subject |
|---|---|
| `cf02836b` | `lane9(v0.4.5): the DISTINCT_donors class is not a soundness class` |
| `c804044f` | `test(lane9): state the proper-rotation guard's scope instead of implying it` |

Files:

| file | change |
|---|---|
| `tools/wrong_donor_groundtruth.py` | **new** (341 lines) — the geometric ground-truth instrument, with the degeneracy argument in its docstring |
| `tools/slot_drift_mechanism.py` | modified (+123/−28 region) — verdict renamed to `distinct_donors_LOCAL` with a warning; fragments paired by body text; haptic donors compared on the whole atom set |
| `tests/unit/test_wrong_donor_verdict.py` | **new** (162 lines at introduction, +28/−13 in `c804044f`) |
| `docs/WRONG_DONOR_v0.4.5.md` | **new** (230 lines) |
| `src/**` | **untouched** |

**Levers: none introduced, none changed.** Levers-OFF bytes unchanged by construction.

Guards run green at land time, as recorded in `cf02836b`: `test_facmer_key` OK (levers off **and**
on), `test_isomer_divergence` OK (off **and** on), `test_regression_stability` goldens OK,
`test_canonical_slots` + `test_canonical_slot_invariance` **OK 46**. Suite baseline at the time:
**674 OK / 3 skip / 3 xfail.** Ruff clean.

Reproduce:

```bash
# ground truth over all 7 (or one molecule); GT_RAW=1 reports the raw fit vertices
PYTHONPATH=src OIN_CANONICAL_BODY=1 OIN_CANONICAL_PERCEPTION=1 \
OIN_CANONICAL_SLOTS=1 OIN_STABLE_METAL_AC=1 GT_RAW=1 \
    .venv/bin/python tools/wrong_donor_groundtruth.py all

# the guard, which needs no dataset
PYTHONPATH=src .venv/bin/python -m unittest tests.unit.test_wrong_donor_verdict -v
```

Cross-check data referenced by the write-up:
`results-v0.4.5-lane2/stacked/canonicality_probe_all.json` (Lane 2's probe output).

---

## Open questions / for the next agent

1. **The §7a wider fold is still the right fix, and it is still a product call.** Lane 9 removed the
   "it would over-fold your 7" objection. The remaining two — *it folds past the geometry's own
   symmetry* and *it needs a Δ/Λ guard* — stand, and the status doc is explicit that the
   recommendation is **not** to take it unilaterally.
2. **Any determinism fix must be aimed at all three tie-break sites, not just the Kabsch one.**
   Fixing `_map_to_template` alone closes at most 3 of the 7 and 0 of the mechanism behind the 23.
3. **Do not resurrect a soundness claim from string evidence.** `distinct_donors_LOCAL` is a local
   heuristic by construction; `test_module_does_not_reintroduce_the_soundness_claim` will fail if
   the old name comes back. If a future pair genuinely looks like a soundness defect, the
   instrument to reach for is `tools/wrong_donor_groundtruth.py`, or a string-level implementation
   of its group-membership test.
4. **The 2 `unparsable` residuals (borane-cluster class) have never been classified at atom level.**
   That is a real, untouched gap in the 32.
5. **PBP's stored template is only approximately C₅-symmetric** (7 decimals; orbit splits at
   ~1e-7). `IMACOO` is decided by that rounding rather than by the molecule, and it is also the one
   molecule Lane 2's ungated rotation-group unification moved. If PBP handling ever matters,
   increase the stored precision or derive the vertices analytically.
6. **Lane 2's "`postpass_BUG_diverges` = 0" is arm-scoped.** It holds for the `on`/`stacked` arms
   only; the levers-OFF and Lane-1-only arms read 4 and 3. Quote it with the arm attached.
7. **The three "resonance pair" molecules are a standing hazard for any symmetry-based classifier.**
   `VEXHIR`/`KISQAG` (acac), `IMACOO` (β-diketiminate) and `HAVGIW` (dipyrrin) have donors that are
   inequivalent **only in the written Kekulé localization**. Any tool that asks
   `CanonicalRankAtoms` whether two donors are equivalent is asking about a localization, not about
   the molecule — expect false "inequivalent" verdicts there.
