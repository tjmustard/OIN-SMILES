# Lane 6 — Metal-locked donor amine (Y1 blind spot **P3**)

**The blind spot:** a secondary amine that is stereogenic *only because a metal occupies its
fourth position* lost its configuration entirely — `POJJOP` and its mirror image encoded
**byte-identically**, with the amine written as a bare `[NH]{0}` carrying no handedness. This is
the user's original motivating case.

---

## ELI5

A nitrogen atom with three different things attached plus a lone pair normally flips inside out
many times a second, like an umbrella in the wind, so it has no fixed handedness and chemistry
software correctly refuses to give it one. But bolt that nitrogen to a metal and the metal takes
the fourth position: the umbrella can no longer turn inside out, and the nitrogen becomes a real
handed centre with two distinct **mirror-image** forms — same atoms, same connections, but like
a left and a right glove, no rotation ever makes one match the other. OIN threw that handedness
away, because the step that splits the complex into pieces removes the metal first, and the
moment the metal is gone the software sees an ordinary flipping amine and erases the label. The
fix is a narrow exemption: keep the label for metal-bonded nitrogens **only**, since the general
erasing rule is right for everything else and protects about a fifth of the collection from
being given handedness it does not have.

---

## The work, visually

```
 WHERE THE CONFIGURATION IS LOST
 ===============================

  input .xyz                metal PRESENT: N has 4 neighbours (M, H, tolyl, CH2-pyridyl)
      │                                    => genuinely stereogenic
      ▼  get_tmc_mol()
      │
      ▼  fragment rebuild  ── metal STRIPPED ──►  N now trivalent (total_degree = 3)
      │                                            │
      │                                            ├─ RDKit: AssignStereochemistry(cleanIt=True)
      │                                            │   clears trivalent N UNCONDITIONALLY
      │                                            │   (correct: a FREE amine really inverts)
      │                                            └─ Zone-A clear in
      │                                                ChiralityRecoveryUtility.recover
      │                                                drops any P/N with total_degree < 4
      ▼
   [Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}     <-- base
   [Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}     <-- mirror.  IDENTICAL.

 THE FIX: A RESTORE AFTER THE LOOP, NOT A NARROWED PREDICATE
 ===========================================================

   perception_tmc.py:1398-1412   plan_locked_donors(tmc_mol)  <- ONCE per molecule, on the
        │                   4-condition gate            metal-PRESENT mol + PRISTINE conformer
        │                                               (before _align_to_pai, which may REFLECT)
        ▼
   perception_tmc.py:1651-1658   stamp_locked_donor_stereo(frag, ...)
        │                   sets an atom PROPERTY (_OIN_LOCKED_STEREO = "CW"/"CCW")
        │                   NOT a tag: properties survive sanitisation; N tags do not
        ▼
   chirality.py:828       restore_locked_donor_tags(rw)   <- LAST action of recover()
        │                   property -> chiral tag, ONLY on atoms left CHI_UNSPECIFIED
        ▼
   inline.py:239          if fragment has "[N@" and lever on:
        │                     mol.SetIntProp("_StereochemDone", 1)
        │                   ^ inline.py reparses each fragment with sanitize=False to attach
        │                     {slot} markers; the following MolToSmiles re-runs
        │                     assignStereochemistry(cleanIt=True) unless _StereochemDone --
        │                     deleting the descriptor at the LAST possible moment
        ▼
   [Pd_SPL].Cc1ccc([N@@H]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}   <-- base
   [Pd_SPL].Cc1ccc([N@H]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}    <-- mirror.  DIVERGED.

 TWO FRAMES, BECAUSE THE @/@@ CHARACTER MEANS TWO DIFFERENT THINGS
 =================================================================
   MEASURED by holding one molecule fixed and feeding the writer two bond orders:

   fragment atom class          does emitted @/@@ depend on bond order?
   ---------------------------  --------------------------------------
   3 or 4 real bonds            YES -- the tag is a PARITY vs bond order
   2 real bonds + implicit H    NO  -- the tag is an ABSOLUTE label

   A metal-locked donor lands in BOTH: PR3 keeps 3 heavy nbrs after the metal
   is cut; a 2-degree amine keeps only 2.

   3+ bonds  ->  parity over first three neighbours IN FRAGMENT BOND ORDER
                 (RDKit's own assignChiralTypesFrom3D rule; the numbering
                  dependence CANCELS -- reorder and both parity and the writer's
                  reading flip together).  Same frame as oin/stable_stereo.py.
   2 bonds   ->  drop the METAL (the neighbour the fragment is missing, standing
                 in for the lone pair) and order the remaining three by
                 CanonicalRankAtoms(breakTies=False) SYMMETRY-CLASS rank.
                 The gate guarantees those three are in distinct classes, so
                 NO index-dependent tiebreak is ever consulted.
                 MEASURED: +9.405 across 6 renumberings, -9.405 for the mirror.

 ⚠ THE HEADLINE PROBLEM: P3 IS NOT USABLE IN THE SHIPPED DEFAULT
 ===============================================================

   OIN_CANONICAL_BODY  (promoted DEFAULT-ON in v0.4.5)
        │
        ▼  canonical_body_emit reparses the ligand body
        ▼  sanitizing a METAL-FREE fragment runs AssignStereochemistry(cleanIt=True)
        ▼  which strips [N@] off a 2-degree amine as a freely inverting nitrogen
        ▼
   the tag is STAMPED and then DISCARDED.
   tests/unit/test_locked_donor.py therefore runs with OIN_CANONICAL_BODY=0.

 THE OBVIOUS FIX WAS TRIED IN v0.4.6 AND IS MEASURABLY WRONG
 ===========================================================
   "copy the chiral tag onto the reparsed donor -- _reparse_once's Guard 2 already
    proves donors[k] <-> new_donors[k] (same element, same heavy degree)"
   4 lines.  P3 emits.  POJJOP PASSES.
        │
        ✗ setting a tag AFTER the sanitize introduces a stereocentre the canonical
          ranker did not account for  ->  moves the canonical WRITE ORDER
          ->  and @/@@ is a PARITY RELATIVE TO THAT ORDER
        │
        ▼  on RIFGUJ_comp_2 the three ring-CARBON tags then flip base vs mirror
        ▼  and the geometry says they MUST NOT:
              AssignStereochemistryFrom3D + rdCIPLabeler give them lowercase `s`
              (PSEUDO-ASYMMETRIC = a RELATIVE all-cis descriptor)
                  atom   6  base=s  mirror=s   same
                  atom   8  base=s  mirror=s   same
                  atom  10  base=s  mirror=s   same
        │
        ▼  REVERTED. Single-centre POJJOP could NOT catch this; the multi-centre
           RIFGUJ guards did. Now pinned by
           test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric

 LEGEND
 ======
   ──► / ▼   flow / sequence in time        ✗  refuted        ⚠  trap or obligation
   {n}       OIN slot marker                @/@@  the SMILES tetrahedral chirality characters
   2-degree  secondary (i.e. N bonded to two carbons plus one H)
   `s`       lowercase CIP label = pseudo-asymmetric (a RELATIVE descriptor)
```

---

## Initial assumptions and hypothesis

`docs/agentic-notes/injectivity/INJECTIVITY_Y1_P3_AMINE.md` had already corrected the *first* hypothesis before the lane
started, and that correction matters: the obvious candidate — the documented `JUCCUH`
`[N@@H]→[NH]` residual — is **not** an encoder blind spot. `convert(JUCCUH)` emits
`...C[N@@H](C)Cc2ccc...` and its mirror emits `[N@H]`, so the encoder **does** distinguish them.
JUCCUH's stereogenic N is a **pendant** backbone amine; its loss happens on the *generated*
re-encode, making it a false-negative / generator effect. **The true P3 blind spot is narrower: a
metal-BOUND secondary amine, stereogenic only because the metal locks its fourth position.**

The lane charter (`spec/handoffs/v0.4.5/Lane6-bound-amine.md`) then made four assumptions, and
**all four were measured wrong**:

| charter assumption | outcome |
|---|---|
| narrow the `total_degree < 4` predicate at `core/chirality.py:722-727` | **would have achieved nothing** — see D1 |
| Lane 2's canonical slot ordering is a dependency | **not a dependency**; no canonical slot permutation is consumed |
| `tests/unit/test_injectivity_probes.py::test_metal_bound_amine_should_diverge_in_raw_string` will flip from xfail to pass | **stays xfail** — that is now a statement about the *default* |
| trivalent P is a real gap with wider blast radius than the amine | **the P gap does not exist as described** — see D5 |
| the motif prevalence is 2.85% (171 / 6000) | **an over-count** — the pre-filter did not test symmetry-distinctness |

---

## What was actually found

### Confirmed

| claim | measurement | source |
|---|---|---|
| P3 is total **encoder** blindness (not merely key blindness) | `convert(POJJOP)` and `convert(mirror)` byte-identical; `raw_equal = True`, `key_equal = True` | `docs/agentic-notes/injectivity/INJECTIVITY_Y1_P3_AMINE.md` |
| the two are genuinely distinct isomers | min proper-rotation mirror RMSD **2.35 Å** | same |
| the configuration IS in the input 3D | signed tetrahedral volume at the locked N: **−9.4** base vs **+9.4** mirror | `docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md` |
| RDKit will not perceive it for us | `AssignAtomChiralTagsFromStructure` returns `CHI_UNSPECIFIED` for a 3-coordinate nitrogen **even with all three neighbours explicit as real atoms with real coordinates** | measured, `locked_donor.py` docstring |
| the standard `[N@]`/`[N@@]` tag IS emittable, if applied after the last sanitise | see the RDKit table below | measured on rdkit **2025.09.3** |
| the descriptor is reproducible in the 2-bond frame | signed volume **+9.405** identical across **6** random renumberings; **−9.405** for the z-mirror | `bac039ef` |
| the over-sensitivity gate is doing real work | `CisPlatin`, `PtMeNH3ClBr` (ammine, M,H,H,H) and `Cis-PtCl2(en)` (primary amine, M,H,H,C) yield **no eligible donor and no string change** | `test_locked_donor.py::TestOverSensitivityGuard` |
| the round-trip key **folds** the descriptor | POJJOP and its mirror still share a key with the lever ON, because the key re-parses each fragment through RDKit | `test_round_trip_key_still_folds_the_descriptor` |

**The RDKit measurement table that drove the notation decision** (rdkit 2025.09.3):

| step | nitrogen | phosphorus (3 distinct substituents) |
|---|---|---|
| `MolToSmiles` of a manually-set tag | **emits** `[N@@H]` | emits `[P@@H]` |
| `AssignStereochemistry(cleanIt=True)` | **clears** | keeps |
| sanitising `MolFromSmiles` | **clears** | keeps |
| `MolFromSmiles(sanitize=False)` | keeps | keeps |
| two identical substituents | n/a | correctly cleared (`CP(C)c1ccccc1`) |

So the standard tag *is* usable, **provided it is applied after the last sanitising step** —
which is why the mechanism stamps an atom **property** during the fragment rebuild and converts
it to a tag as `recover()`'s final action. **Properties survive sanitisation; tags on nitrogen do
not.**

### Measured prevalence and blast radius

400 corpus molecules, seed 42, `cat` + `photo`, perception + encode only (no 3D generation, so
no timeout confound), 25 s per-molecule cap. 395 encoded, 2 errors, 3 timeouts.

| | count | share of 395 |
|---|---:|---:|
| at least one eligible metal-locked donor | **21** | **5.3%** |
| — with an eligible **N** | 18 | 4.6% |
| — with an eligible **P** | 3 | 0.8% |
| **OIN string changes with the lever ON** | **18** | **4.6%** |
| eligible N centres / eligible P centres | 27 / 7 | |

All 18 changed molecules have `nP = 0`, and no molecule had both an eligible N and an eligible P
(18 + 3 = 21), so the three eligible-but-unchanged molecules are exactly the three P-bearing
ones. Hence the single most useful number in the lane:

> **All 7 eligible metal-locked phosphorus donors, across all 3 molecules carrying them, already
> had a chiral tag from the existing Zone-A lone-pair path. The restore was a no-op for every
> one.**

### Encoder results, lever ON (`mirror` = z-reflection)

| molecule | lever OFF | lever ON |
|---|---|---|
| `POJJOP` (1 locked amine, sole stereocentre) | mirror byte-identical | `[N@@H]{0}` vs `[N@H]{0}`; **`mirror == swap(base)`** |
| `RIFGUJ_comp_2` (3 Cu-bound amines) | mirror byte-identical | all three invert, **nothing else moves** |
| `TEGFET`, `KANYUW`, `LARYEI` | already diverged via a second stereocentre | still diverge, now with explicit amine descriptors |
| `JUCCUH` | pendant `[N@@H]` already captured | **byte-identical** — no eligible donor |

### Refuted — see *Dead ends*. Five distinct refutations, four of them of the lane's own charter.

---

## What was done

Five commits on `swimlane/v045-lane6`, merged into `release/v0.4.5` as `df7417a9`:

| commit | title |
|---|---|
| `bac039ef` | `feat(lane6): recover metal-locked N/P donor chirality behind OIN_EMIT_LOCKED_DONOR` |
| `b6e2310a` | `test(lane6): three-property + over-sensitivity guards for the locked-donor descriptor` |
| `9bc8b025` | `docs(lane6): the P3 write-up, plus the key-fold guard and the xfail correction` |
| `d64a4154` | `docs(lane6): measured prevalence, and the trivalent-P gap refuted` |
| `63ba9090` | `docs(lane6): record the measured acceptance table` |
| (`fc2265ee`) | `docs(v0.4.5): Lane 6 closes P3 -- and corrects three things I asserted` |

### The carve-out: metal-bonded N (and P) ONLY

`src/oinsmiles/oin/locked_donor.py` is new; the touched existing files are
`src/oinsmiles/core/chirality.py`, `src/oinsmiles/utils/perception_tmc.py` and
`src/oinsmiles/oin/inline.py`.

**The general clearing rule is left UNCONDITIONAL on purpose**, and this is the load-bearing
design decision. The unconditional clear is what keeps a genuinely invertible free amine and a
symmetric phosphine from acquiring spurious sp3 handedness (see the `with_stereo=False` note at
`core/translator.py:105-116`), and **that protection covers roughly a fifth of the corpus.
Widening the exemption beyond metal-bound N will regress them.**

So the carve-out is a **restore that runs after the whole existing loop**
(`core/chirality.py:828`, the last stereo action of `ChiralityRecoveryUtility.recover`), reading
a property stamped during the fragment rebuild. It writes **only** atoms whose tag is currently
`CHI_UNSPECIFIED`, so the Zone-A lone-pair P branch and the ≥4-neighbour verify-and-flip both
keep priority. Blast radius is exactly *"metal-locked donors that would otherwise carry no tag at
all"* — which is also precisely the gap Lane 8 could not reach, since
`stable_stereo.py:112-118` only corrects tags that **already exist**.

### The stereogenicity gate — `plan_locked_donors()`

Four conditions, all narrow on purpose:

1. the atom is **N or P** and is **not aromatic** (an aromatic nitrogen cannot be a tetrahedral
   stereocentre, and a pyridine donor must never be touched);
2. bonded to **exactly one** transition metal (a bridging donor's configuration is not a property
   of one centre);
3. **exactly four** neighbours in the metal-present mol — the metal plus three, which is what
   "the metal locks the fourth position" *means*;
4. those four in **four distinct symmetry classes** via `CanonicalRankAtoms(breakTies=False)`.

**Condition 4 is the over-sensitivity guard**, so a metal-bound **ammine** (M,H,H,H) or a
**primary** amine (M,H,H,R) emits **nothing**. *That guard exists because the axial lane had to
fix exactly that failure once already* — an encoder claiming stereochemistry that does not exist
is a defect of the same magnitude as the blind spot it was meant to fix.

`breakTies=False` is not an implementation detail: `breakTies=True` was measured **not**
invariant under input renumbering (**2–11 distinct rank vectors over 20 renumberings**), and that
is pinned by `TestPlanIsRenumberingInvariant` so nobody "improves" it later.

The gate also correctly declines **three corpus molecules the Y1 write-up had listed as blind**:

| molecule | why it is not a stereocentre |
|---|---|
| `CULGOF` | Rh-bound morpholine N — the two ring α-carbons are symmetry-equivalent |
| `ATEFIP` | Pd-bound diethylamine N (two equivalent ethyls) and a primary `[NH2]` donor |
| `MOKCEV` | Ir-bound macrocyclic triamine N with equivalent ring carbons |

Also `_PLANARITY_TOL = 0.05` on the |triple product| of the three unit reference vectors: below
that the centre is too close to planar to read a handedness from, so **nothing is stamped rather
than resolving it by numerical noise**. Same value and rationale as
`oin/stable_stereo.py::_PLANARITY_TOL`.

### Absolute sense is an OIN convention, stated as such

RDKit has no answer to match — `AssignAtomChiralTagsFromStructure` returns `CHI_UNSPECIFIED` for
a 3-coordinate nitrogen even with all three neighbours explicit. So "which of `@`/`@@` is R" is
an **OIN convention** fixed by the rule above, not an inherited RDKit or CIP one. That is enough
for injectivity — reproducible, orientation-invariant, and inverting under reflection — and all
three are pinned by tests so the convention cannot drift silently.

### Why the mirror assertions compare whole strings, never tag counts

This is a methodological decision worth carrying forward. **Counting is blind to a symmetric
swap:** three `@@` plus three `@` mirrors to three `@` plus three `@@` — *identical counts*. An
earlier hand-run of exactly this check reported a **false pass** for that reason.

So `test_locked_donor.py` compares **whole strings under an `@`↔`@@` swap**. And it needs two
swap functions, because the strict form is wrong in general:

- `_swap_tags(s)` — exchange *every* `@`/`@@`. Used on **POJJOP**, whose whole string is
  byte-stable under renumbering with the lever off.
- `_swap_n_tags(s)` — exchange only `[N@…]`↔`[N@@…]`. Used on **RIFGUJ_comp_2**, because
  `mirror == swap(EVERY tag)` is **false** there, correctly: its 1,3,5-trisubstituted cyclohexane
  **ring-carbon tags encode a RELATIVE (all-cis) arrangement and are measured NOT to invert under
  reflection** — verified with the lever **off**, so it is not something this lane introduced.
  (Lane 8 recorded the same shape of problem for an η-arene whose mirror differed only in the
  winding character.) The assertion is therefore **"every locked-nitrogen descriptor inverted and
  nothing else in the string moved"** — still a whole-string comparison, still immune to
  count-blindness.

**Full-string renumbering stability cannot be asserted on RIFGUJ.** With the lever *off* it
drifts in **3 of 3** renumberings — ring-carbon tags flip, one loses its tag entirely, and the
`{2}`/`{3}` slot numbers swap. That is the pre-existing 13% stereo-flip class (Lane 8's
`OIN_STABLE_STEREO`) plus slot drift (Lane 2). Asserting it here would import two other lanes'
open defects into this lane's guard, so the renumbering assertion is scoped to the nitrogen
descriptors while rotation and reflection — both clean at whole-string level — keep the strict
form.

### Why an out-of-band token was rejected

An out-of-band sidecar (`|amine:0+|`, following the landed `|ax:±|`) **would** survive an RDKit
re-parse, which is the one real advantage. Rejected anyway, for two reasons and not for effort:

1. it needs its own atom identity, sign convention and multi-token ordering — **three fresh
   chances to re-run the Y2 "sorted by sign, therefore reflection-invariant" mistake**;
2. it would give phosphorus a **second** representation on top of the `[P@]` the existing Zone-A
   lone-pair path already emits.

---

## Dead ends and refutations

### D1 — "Narrow the `total_degree < 4` predicate at `core/chirality.py:722-727`"

This is what the charter prescribed. **Killed by reading the two lines above it.** `recover()`
calls `Chem.AssignStereochemistry(cleanIt=True, force=True)` *before* the Zone-A branch, and that
call clears a trivalent nitrogen **regardless** of what the branch then does. **Exempting the
atom from the Zone-A clear exempts it from nothing.**

And the unconditional clear must stay unconditional independently of that: it is what protects
the ~20% of complexes that would otherwise acquire spurious sp3 handedness.

**Replacement:** a restore *after* the whole loop, writing only atoms left `CHI_UNSPECIFIED`.
The comment at `core/chirality.py:795-803` now records this so the "obvious" narrowing cannot be
re-attempted.

### D2 — "Use `bound_amine_centers`' neighbour ordering"

`tools/injectivity/config_oracle.py::bound_amine_centers` orders the donor's neighbours by
`CanonicalRankAtoms(breakTies=True)`.

**Killed by:** that ordering is **not invariant under input renumbering** — measured **2–11
distinct rank vectors over 20 renumberings**. Fine as an oracle on one fixed molecule; **useless
as an emission mechanism.** Not used here.

### D3 — "Use the fragment's own bond order and let the canonical writer normalise it"

This is Lane 8's argument for `restamp_fragment_chirality`, it is **correct there**, and it was
the lane's **first implementation**.

**Killed by:** it drifted under renumbering in **2 of 3 trials on POJJOP**. The reason, measured
by holding one molecule fixed and handing the SMILES writer two different bond orders at the
stereocentre:

| fragment-atom class | does the emitted `@`/`@@` depend on bond order? |
|---|---|
| 3 or 4 real bonds | **yes** — the tag is a parity relative to bond order |
| 2 real bonds + implicit H | **no** — the same tag emits the same character either way |

A metal-locked donor lands in **both** classes (a PR₃ phosphine keeps three heavy neighbours once
the metal bond is cut; a secondary amine or phosphine keeps only two), so **two frames** are
needed. The fix is described under *What was done*; the 2-bond frame's key property is that the
stereogenicity gate has already guaranteed the three reference neighbours sit in **distinct**
symmetry classes, so **no index-dependent tiebreak is ever consulted** — which is exactly what
separates it from D2.

### D4 — "The xfail will flip to passing"

`tests/unit/test_injectivity_probes.py::test_metal_bound_amine_should_diverge_in_raw_string`
remains an `@expectedFailure`.

**Not killed by a bug — killed by the acceptance criterion.** The fix is gated, and
"levers-OFF output is byte-identical" is the *harder* requirement. The xfail is now a statement
about the **default**, and its docstring says when to flip it (when the lever is promoted). The
capability is guarded permanently in `tests/unit/test_locked_donor.py` instead.

### D5 — "Trivalent P is a real gap"

`docs/agentic-notes/v0.4.5/RENUMBERING_INSTABILITY_v0.4.5.md` handed Lane 6 the trivalent-P case on the reasoning that
`stable_stereo.py:112-118` can only correct tags that already exist, so a P donor cleared by the
Zone-A rule has nothing to restamp. **The reasoning is sound; the case is not there.**

**Killed by:** the 400-molecule prevalence scan (all 7 eligible P donors already carried a tag),
plus direct measurement on the two molecules named as Lane 8's residuals:

- **`FEQFIS_comp_0`** — its metal-bound phosphorus is P(N)(O)(O)Au, **aromatic**, and its four
  neighbours fall in only **3 distinct symmetry classes** (the two oxygens are equivalent). It is
  **not** a stereocentre; emitting for it would be the over-sensitivity failure, not a fix. Its
  renumbering drift is in a `[C@@H]` carbon anyway, **not at P**.
- **`CEBVIR_comp_0`** — four aromatic nitrogen donors, each with **3** neighbours
  (pyridine-type), so none is a tetrahedral stereocentre either.

Both correctly declined; the lever changes neither string. **The P half of this lane is
implemented, tested, and measured to be unnecessary on this corpus. Do not quote this lane as
having closed a phosphorus defect.**

### D6 — "The motif prevalence is 2.85%"

Y1's figure (**171 / 6000**) came from a pure-geometry pre-filter — "N within 2.6 Å of a
transition metal with exactly 1 H and 2 C neighbours" — which **does not test
symmetry-distinctness**. **Killed by** the gate declining CULGOF, ATEFIP and MOKCEV, three
molecules the Y1 write-up had listed as blind. So 2.85% is an **over-count of the motif that is
actually stereogenic**. The measured replacement is 5.3% eligible / 4.6% string-changing on a
different (400-molecule) sample — note the two numbers are **not directly comparable**, being
different populations and different predicates.

### D7 — ⚠ **THE BIG ONE: "copy the chiral tag onto the reparsed donor"** (tried in v0.4.6, MEASURED WRONG)

**What it was.** Four lines. The correspondence is already available:
`oin/canonical_body.py::_reparse_once`'s **Guard 2** already proves `donors[k] ↔ new_donors[k]`
(same element, same heavy degree), so copying the chiral tag from `mol` onto `reparsed` looks
safe. **It makes P3 emit under `OIN_CANONICAL_BODY`, and POJJOP passes.**

**Why it is still wrong.** Setting a chiral tag **AFTER** the sanitize introduces a stereocentre
**the canonical ranker did not account for**, which moves the canonical **WRITE ORDER** — and
`@`/`@@` is a **parity relative to that order**, not an absolute label.

**The measurement that killed it.** On `RIFGUJ_comp_2` (three Cu-bound amines on one cyclohexane)
the three ring-**CARBON** tags then flip between a structure and its mirror. The geometry says
they must not: `AssignStereochemistryFrom3D` + `rdCIPLabeler` label those carbons lowercase
**`s`** — **pseudo-asymmetric**, i.e. a **RELATIVE** (all-cis) descriptor — and they read `s`
identically for the structure **and** its reflection:

```
atom   6  base=s  mirror=s   same
atom   8  base=s  mirror=s   same
atom  10  base=s  mirror=s   same
```

So the restoration **silently rewrote stereochemistry that must not move.**

**Who caught it.** `TestMultiCentreDescriptor::test_flips_under_reflection` and
`TestLeverOnDivergesOnEnantiomers::test_three_locked_amines_all_invert_together`. **Single-centre
POJJOP could not** — the Y2 lesson intact: a guard that only exercises the easy case confirms
wrong beliefs.

**Disposition.** Reverted. The mechanism is written into three places so it cannot be silently
re-attempted: `oin/canonical_body.py::_reparse_once` (lines ~245-259, an explicit "DO NOT restore
… Tried in v0.4.6 and MEASURED WRONG" block), `oin/levers.py::_HELD_OFF["OIN_EMIT_LOCKED_DONOR"]`,
and `tests/unit/test_locked_donor.py::_LeverBase`'s docstring. Plus a new **lever-independent**
guard, `test_locked_donor.py::TestRifgujRingCarbonsArePseudoAsymmetric`, which runs in the
**shipped configuration** (all v0.4.5 defaults, locked-donor lever off) so it cannot be pinned
away.

**A correct fix must preserve the tag WITHOUT perturbing the ranking.** Two named candidates:
keep the donor **bracketed** through the sanitize, or **re-derive parity from the parent
geometry** once the write order is fixed.

**Also worth noting:** v0.4.5's decision to *defer* this rather than rush it was therefore
correct.

### D8 — Known residual, honestly unresolved: internally compensated donor pairs

`ABIFAV_comp_0` (Pd, two adamantanecarboxylate O donors and two *N*-isopropylcyclohexylamine
donors) emits **both** descriptors — `{1}` as `[N@@H]` and `{3}` as `[N@H]` — and yet base and
mirror remain **byte-identical** with the lever ON. The two amine ligands are constitutionally
identical and carry **opposite** configurations, so the pair is **internally compensated**:
reflecting the molecule exchanges which ligand is which, and the descriptor pair is
reflection-invariant. On that reading the identical encoding is *correct*, and the rigid-mirror
oracle's "distinct, RMSD 3.30 Å" verdict is the known conformation-inflated artefact for flexible
ligands. `LARYEI_comp_0` shows the same pattern at its two ethylenediamine-derived nitrogens, and
there the string still differs because its backbone carbons do invert.

**What is NOT established:** whether the fold in `ABIFAV` is purely the (R,S) compensation
described above, or whether the reflection-sensitivity of **slot assignment** is also
contributing. Distinguishing them needs a **conformer-independent isomer oracle**, which this
lane did not build. Recorded rather than resolved.

---

## Where it landed

**Lever:** `OIN_EMIT_LOCKED_DONOR` — **default OFF** (in `levers.py::_HELD_OFF`). Env var name in
code as `locked_donor.ENV_LEVER`. Routed through `levers.lever_enabled`, so
`OIN_EMIT_LOCKED_DONOR=0` **disables** it; the bare `bool(os.environ.get(...))` this replaced did
the opposite, because `"0"` is a non-empty string.

**Code:**

| stage | location |
|---|---|
| eligibility + descriptor | `src/oinsmiles/oin/locked_donor.py` — `plan_locked_donors`, `_reference_neighbours`, `_tag_name_from_geometry`, `stamp_locked_donor_stereo`, `restore_locked_donor_tags` |
| plan computed once, on the **pristine** conformer | `src/oinsmiles/utils/perception_tmc.py:1398-1412` (before `_align_to_pai`, which can reflect the coordinates and would invert the recovered sign) |
| property stamped during fragment rebuild | `src/oinsmiles/utils/perception_tmc.py:1651-1658`, carrying `LOCKED_TAG_PROP = "_OIN_LOCKED_STEREO"` (`"CW"`/`"CCW"`) |
| property → tag, as `recover()`'s LAST action | `src/oinsmiles/core/chirality.py:828` |
| survive the slot-marker round trip | `src/oinsmiles/oin/inline.py:239` — `_StereochemDone` marked, **gated twice** (lever on **AND** the fragment actually contains `[N@`) so no other fragment's serialization changes |

**Guard tests — `tests/unit/test_locked_donor.py`, 16 tests, all OK** (measured on `main`):

| class | tests | what it pins |
|---|---:|---|
| `TestDefaultOff` | 2 | POJJOP mirror byte-identical with the lever off; bare `[NH]{0}`; no `[N@` on RIFGUJ |
| `TestLeverOnDivergesOnEnantiomers` | 3 | POJJOP enantiomers diverge with `mirror == _swap_tags(base)`; **`test_three_locked_amines_all_invert_together`**; **`test_round_trip_key_still_folds_the_descriptor`** |
| `TestDescriptorIsCanonical` | 3 | three-property test on POJJOP at **whole-string** strictness — 3 renumberings, 3 proper rotations, reflection |
| `TestMultiCentreDescriptor` | 3 | the same on RIFGUJ, scoped per the reasons above — descriptors over 4 renumberings; whole-string rotation; whole-string reflection under `_swap_n_tags` |
| `TestOverSensitivityGuard` | 3 | ammine / primary amine not eligible **and** no string change; BINAP as a P-donor control; POJJOP as the positive control so an over-tight gate cannot pass silently |
| `TestPlanIsRenumberingInvariant` | 1 | eligibility itself is renumbering-invariant (RIFGUJ: 3 donors on 4 of 4 renumberings) |
| `TestRifgujRingCarbonsArePseudoAsymmetric` | 1 | **the D7 guard.** Runs in the SHIPPED configuration, lever-independent |

Note the arithmetic: the lane's own acceptance table recorded **620 OK / 3 skip / 3 xfail**
against a **605 / 3 / 3** baseline on the same branch — delta **+15**, exactly the new file, so
every delta accounted for and zero regressions. The file now has **16** tests; the 16th is
`TestRifgujRingCarbonsArePseudoAsymmetric`, added in v0.4.6 when D7 was refuted.

**Fixtures:** `tests/fixtures/POJJOP.xyz` (single Pd-bound 2° amine — the motivating case) and
`tests/fixtures/RIFGUJ_comp_2.xyz` (**three** Cu-bound 2° amines on one cyclohexane — the
multi-centre case, and the reason the mirror check compares whole strings). Controls:
`CisPlatin.xyz`, `PtMeNH3ClBr-Cis.xyz`, `Cis-PtCl2(en).xyz`, `PdCl2-R-BINAP.xyz`.

### ⚠ What is NOT usable in the shipped default, and exactly why

**P3 is NOT usable in the shipped default configuration.** Two independent reasons:

1. **`OIN_EMIT_LOCKED_DONOR` is default OFF** — the standard information-ADDING trade. With the
   lever unset nothing is stamped, so the restore loop in `recover()` and the `_StereochemDone`
   guard in `inline.py` are both no-ops over an empty set, and output is byte-identical.
2. **Even with the lever ON, it is INCOMPATIBLE with `OIN_CANONICAL_BODY`, which was promoted
   default-ON in v0.4.5.** `canonical_body_emit` reparses the ligand body, and sanitizing a
   **metal-free** fragment runs `AssignStereochemistry(cleanIt=True)`, which strips the chiral tag
   off a 2° amine as a freely inverting nitrogen — **the exact RDKit behaviour this descriptor
   exists to work around.** With both levers on, the tag is stamped upstream and then **discarded**.
   `tests/unit/test_locked_donor.py::_LeverBase` therefore pins `OIN_CANONICAL_BODY=0`, and says
   so in its docstring — *stated rather than hidden, because the pin makes this suite green
   against a configuration nobody ships.*

**A third, separate limitation even in the lever-ON + canonical-body-OFF configuration:** for
nitrogen the descriptor is **encoder-authoritative but not RDKit-round-trippable**. A sanitising
re-parse of the OIN drops it, so the generator cannot rebuild it and the round trip reports a
**loud false negative** rather than a silent collision — the same trade the axial lever makes.
And the round-trip comparison key **folds** the descriptor (measured: POJJOP and its mirror share
a key with the lever ON), which cuts both ways and both ways are load-bearing:

- promoting the lever **cannot move** `facmer_divergent` or any other harness count;
- the harness **cannot confirm** the descriptor either. **Only the raw string can.**

**Generator support was not attempted**, and that is the concrete blocker to promoting the lever.
RDKit drops `[N@]` on the parse the generator performs, so selection-by-descriptor (the
`_axial_narrow` precedent) would need the descriptor carried **out of band** into the parsed
structure first. That is a generator-side task.

**Reproduce:**

```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES
export PYTHONPATH=$PWD/src:$PWD
V=$PWD/.venv/bin/python

# the guards
$V -m unittest tests.unit.test_locked_donor

# the defect and the fix, side by side
$V -m tools.injectivity.twin_collision tests/fixtures/POJJOP.xyz
OIN_EMIT_LOCKED_DONOR=1 $V -m tools.injectivity.twin_collision tests/fixtures/POJJOP.xyz
```

### ⚠ Source-document discrepancies

1. **Commit `fc2265ee`'s message names the lever `OIN_EMIT_BOUND_AMINE`.** No such lever exists;
   the shipped name is **`OIN_EMIT_LOCKED_DONOR`**. Treat the commit message as stale on that one
   point.
2. **`docs/agentic-notes/v0.4.5/LOCKED_DONOR_v0.4.5.md` §9 records the test file as +15 tests.** It now has **16** —
   the extra one is v0.4.6's D7 guard. Not a contradiction, just a version skew.
3. **`docs/KNOWN_LIMITATIONS.md`'s P3 bullet still reads "recoverable, deferred to v0.4.5"** and
   describes the emit as requiring "a Zone-A carve-out + canonical ordering + generator". The
   Zone-A carve-out landed (as a restore, not a narrowing), canonical ordering turned out **not**
   to be needed (D3's two-frames finding), and the generator half is still open. The bullet is not
   wrong about the *state*, only about the *shape of the remaining work*.
4. **`docs/agentic-notes/injectivity/INJECTIVITY_Y1_P3_AMINE.md` still cites the 2.85% figure** without the over-count
   correction (D6). `docs/agentic-notes/v0.4.5/LOCKED_DONOR_v0.4.5.md` §5 and §5b carry the correction.

---

## Open questions / for the next agent

### The promotion gate for `OIN_EMIT_LOCKED_DONOR` — in dependency order

1. **Fix the `OIN_CANONICAL_BODY` incompatibility FIRST.** Nothing else matters until the
   descriptor survives the default configuration. **Do not re-try the four-line tag copy** — it is
   refuted (D7) and guarded. The two candidate approaches named in the code are:
   - **keep the donor bracketed through the sanitize**, so the canonical ranker sees the
     stereocentre when it computes the write order rather than after; or
   - **re-derive parity from the parent geometry once the write order is fixed** — i.e. compute
     `@`/`@@` *against the post-reparse canonical order* rather than transplanting a parity that
     was computed against a different one.

   The acceptance test for any attempt is already written:
   `TestRifgujRingCarbonsArePseudoAsymmetric` must stay green **and**
   `TestMultiCentreDescriptor::test_flips_under_reflection` must hold with `OIN_CANONICAL_BODY=1`.
   *Any* candidate fix must be validated on **RIFGUJ**, never on POJJOP alone — POJJOP is
   single-centre and structurally cannot detect a write-order perturbation.

2. **Then the generator.** RDKit drops `[N@]` on the generator's parse, so
   selection-by-descriptor needs the configuration carried out of band into `ParsedOIN` first
   (`generation/oin_parser.py` already strips sidecars, and `original_oin` is already the channel
   the axial and `|mc:|` tokens use). Only then can an `_axial_narrow`-style pass exist. Until
   then, promoting the lever converts a silent false positive into a loud false negative with no
   path to recovery.

3. **Then un-fold the key.** The round-trip key folds the descriptor. Once the generator can
   reproduce it, the fold must be removed in the same commit that promotes the lever — the same
   obligation recorded for `_AXIAL_TOKEN_RE` and `_METAL_CONFIG_TOKEN_RE`. **A key that folds an
   axis is not a valid acceptance predicate for that axis.**

4. **Then flip the xfail.**
   `tests/unit/test_injectivity_probes.py::test_metal_bound_amine_should_diverge_in_raw_string`
   — its docstring already says when.

### Specific open technical questions

- **Is `ABIFAV_comp_0`'s fold pure (R,S) compensation, or is slot assignment contributing?**
  (D8.) This needs a **conformer-independent isomer oracle** — the rigid-mirror oracle
  over-reports for flexible ligands, so its "distinct, RMSD 3.30 Å" verdict cannot settle it. Lane
  7's `tools/injectivity/torsion_oracle.py` is the closest existing instrument and was not applied
  to this case.
- **Do the two prevalence numbers reconcile?** Y1's 2.85% (171/6000, geometry pre-filter, no
  symmetry test) and Lane 6's 5.3% eligible / 4.6% string-changing (400 molecules, full gate) are
  different populations measured by different predicates. Running the *gate* over the *6000*
  sample would give a single comparable rate, and would also tell you the real blast radius of a
  promotion. **Load-independent** (perception + encode only), so it can be taken while a sweep
  runs.
- **The 2-bond frame's sign convention is OIN-local and untied to CIP.** It is reproducible,
  orientation-invariant and reflection-inverting — sufficient for injectivity — but nothing maps
  `@`/`@@` to R/S for a metal-locked nitrogen. If a downstream consumer ever needs a CIP label
  there, that mapping has to be defined and pinned; do not let it be inferred.
- **`_PLANARITY_TOL = 0.05` has no measured corpus distribution.** It is inherited from
  `stable_stereo.py`. Nobody has checked how many eligible metal-locked donors fall near it (and
  therefore silently emit nothing).
