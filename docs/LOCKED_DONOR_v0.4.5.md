# Metal-locked donor chirality — v0.4.5 Lane 6 (injectivity blind spot P3)

> **Headline:** the user's motivating case is closed at the encoder. `POJJOP` — square-planar
> Pd whose sole stereocentre is a Pd-bound secondary amine — and its mirror encoded
> **byte-identically**; behind `OIN_EMIT_LOCKED_DONOR` they now encode as `[N@@H]{0}` and
> `[N@H]{0}`, and the mirror string is the base string with every locked-nitrogen descriptor
> inverted and nothing else changed. Lever **default OFF**; levers-OFF output is byte-identical.

Lane: `swimlane/v045-lane6`. Code: `src/oinsmiles/oin/locked_donor.py`,
`src/oinsmiles/core/chirality.py`, `src/oinsmiles/utils/xyz2mol.py`,
`src/oinsmiles/oin/inline.py`. Tests: `tests/unit/test_locked_donor.py`.

---

## 1. The defect

A secondary amine that is stereogenic **only because the metal occupies its fourth position**
lost its configuration entirely:

```
base   [Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}
mirror [Pd_SPL].Cc1ccc([NH]{0}Cc2ccccn{1}2)cc1.[Cl]{2}.[Cl]{3}   ← identical
```

Two rules cause it, and both are correct in general:

* the Zone-A clear in `ChiralityRecoveryUtility.recover` (`core/chirality.py`) drops the chiral
  tag of any P/N whose `total_degree < 4` once the metal is stripped; and
* RDKit clears trivalent-nitrogen chirality unconditionally, because a free amine really does
  invert.

Neither is true for a metal-locked donor. The metal bond is what makes the centre stereogenic
*and* what stops the inversion.

## 2. Where the carve-out went, and why not where the plan said

The plan said to narrow the `total_degree < 4` predicate at `core/chirality.py:722-727`. **That
would not have worked, and it is left unconditional on purpose.** Two reasons:

1. Two lines earlier `recover()` calls `Chem.AssignStereochemistry(cleanIt=True, force=True)`,
   which clears a trivalent nitrogen regardless of what the Zone-A branch then does. Exempting
   the atom from the Zone-A clear exempts it from nothing.
2. The unconditional clear is load-bearing. It is what keeps a genuinely invertible amine and a
   symmetric phosphine from acquiring spurious sp3 handedness (see the `with_stereo=False` note
   at `core/translator.py:105-116`), and that protection covers a large slice of the corpus.

So the carve-out is a **restore that runs after the whole existing loop**, reading a property
stamped during the fragment rebuild. It writes only atoms whose tag is `CHI_UNSPECIFIED`, so the
Zone-A lone-pair P branch and the `>=4`-neighbour verify-and-flip both keep priority. Blast
radius is exactly "metal-locked donors that would otherwise carry no tag at all" — which is also
precisely the gap Lane 8 could not reach, since `stable_stereo.py:112-118` only corrects tags
that already exist.

## 3. Notation: `[N@]`/`[N@@]` and `[P@]`/`[P@@]`, one notation

Measured on rdkit 2025.09.3:

| step | nitrogen | phosphorus (3 distinct substituents) |
|---|---|---|
| `MolToSmiles` of a manually-set tag | **emits** `[N@@H]` | emits `[P@@H]` |
| `AssignStereochemistry(cleanIt=True)` | **clears** | keeps |
| sanitising `MolFromSmiles` | **clears** | keeps |
| `MolFromSmiles(sanitize=False)` | keeps | keeps |
| two identical substituents | n/a | correctly cleared (`CP(C)c1ccccc1`) |

So the standard tag *is* usable, provided it is applied after the last sanitising step — which is
why the mechanism stamps an atom **property** during the fragment rebuild and converts it to a
tag as `recover()`'s final action. Properties survive sanitisation; tags on nitrogen do not.

One further place had to be taught this. `oin/inline.py` re-parses each fragment SMILES with
`sanitize=False` to attach the `{slot}` markers, and the `MolToSmiles` that follows re-runs
`assignStereochemistry(cleanIt=True)` whenever the mol carries no `_StereochemDone` — deleting
the descriptor at the last possible moment. Marking perception done there fixes it, gated twice
(lever on **and** the fragment actually contains `[N@`) so no other fragment's serialization
changes.

**The honest limitation.** For nitrogen the descriptor is *encoder-authoritative but not
RDKit-round-trippable*: a sanitising re-parse of the OIN drops it. The generator therefore cannot
rebuild it, so the round trip reports a loud false negative rather than a silent collision — the
same trade the axial lever makes, and the reason this lever is default-OFF. Measured: the
round-trip comparison key **folds** the descriptor (`POJJOP` and its mirror still share a key
with the lever ON), so promoting the lever cannot move `facmer_divergent` or any other harness
count — and equally, the harness cannot *confirm* the descriptor. Only the raw string can.

An out-of-band token (`|amine:0+|`, following `|ax:±|`) would survive re-parse. Rejected: it
needs its own atom identity, sign convention and multi-token ordering — three fresh chances to
re-run the Y2 "sorted by sign, therefore reflection-invariant" mistake — and it would give
phosphorus a *second* representation on top of the `[P@]` the existing Zone-A lone-pair path
already emits.

## 4. The reproducible-ordering problem, and the measurement that reframed it

This was the crux, and both of the obvious answers are wrong.

**Wrong answer 1 — `bound_amine_centers`' ordering.** `tools/injectivity/config_oracle.py`
orders the donor's neighbours by `CanonicalRankAtoms(breakTies=True)`, which was measured **not**
invariant under input renumbering (2–11 distinct rank vectors over 20 renumberings). Fine as an
oracle on one fixed molecule; useless as an emission mechanism. Not used here.

**Wrong answer 2 — "use the fragment's own bond order and let the canonical writer normalise
it".** This is Lane 8's argument for `restamp_fragment_chirality`, it is correct there, and it
was my first implementation. It **drifted under renumbering in 2 of 3 trials on POJJOP.**

The reason, measured by holding one molecule fixed and handing the SMILES writer two different
bond orders at the stereocentre:

| fragment-atom class | does the emitted `@`/`@@` depend on bond order? |
|---|---|
| 3 or 4 real bonds | **yes** — the tag is a parity relative to bond order |
| 2 real bonds + implicit H | **no** — the same tag emits the same character either way |

A metal-locked donor lands in *both* classes: a PR₃ phosphine keeps three heavy neighbours once
the metal bond is cut, a secondary amine or phosphine keeps only two. So there are two frames:

* **3+ real bonds** — parity over the first three neighbours in fragment bond order, which is
  RDKit's own `assignChiralTypesFrom3D` rule. The numbering-dependence cancels: reorder the bonds
  and both the parity and the writer's reading of it flip together. Same frame as
  `oin/stable_stereo.py`.
* **2 real bonds + implicit H** — the tag is an absolute label, so the frame must not depend on
  numbering at all. Drop the metal (the neighbour the fragment is missing, standing in for the
  lone pair, i.e. the one RDKit's rule would place opposite the other three) and order the
  remaining three by `CanonicalRankAtoms(breakTies=False)` symmetry-class rank. The
  stereogenicity gate has already guaranteed those three are in distinct classes, so **no
  index-dependent tiebreak is ever consulted** — which is exactly what separates this from wrong
  answer 1. Measured on POJJOP: signed volume **+9.405 across 6 random renumberings, −9.405 for
  the z-mirror**.

Both frames then use the same formula and sign convention as RDKit (`t = (u₀ × u₁) · u₂` over
unit vectors, positive → `CCW`). Only the source of the ordering differs.

**Absolute sense is an OIN convention here, not an inherited one.** RDKit has no answer to match:
`AssignAtomChiralTagsFromStructure` returns `CHI_UNSPECIFIED` for a 3-coordinate nitrogen even
with all three neighbours explicit as real atoms with real coordinates (measured) — which is the
whole reason this code has to exist. Reproducible + orientation-invariant + inverting under
reflection is what injectivity needs, and all three are pinned by tests so the convention cannot
drift silently.

**Lane 2 was therefore not a dependency.** No canonical slot permutation is consumed.

## 5. Stereogenicity gate

Eligibility (`plan_locked_donors`), all four conditions narrow on purpose:

1. N or P, **not aromatic** (an aromatic nitrogen cannot be a tetrahedral stereocentre, and a
   pyridine donor must never be touched);
2. bonded to **exactly one** transition metal (a bridging donor's configuration is not a property
   of one centre);
3. **exactly four** neighbours in the metal-present mol — the metal plus three, which is what
   "the metal locks the fourth position" means;
4. those four in **four distinct symmetry classes** (`CanonicalRankAtoms(breakTies=False)`).

Condition 4 is the over-sensitivity guard, and it is doing real work. Verified emitting nothing:
`CisPlatin` and `PtMeNH3ClBr` (ammine, M,H,H,H), `Cis-PtCl2(en)` (primary amine, M,H,H,C).
Verified from the corpus, where the same gate correctly declines three molecules the Y1 write-up
had listed as blind:

| molecule | why not a stereocentre |
|---|---|
| `CULGOF` | Rh-bound morpholine N — the two ring α-carbons are symmetry-equivalent |
| `ATEFIP` | Pd-bound diethylamine N (two equivalent ethyls) and a primary `[NH2]` donor |
| `MOKCEV` | Ir-bound macrocyclic triamine N with equivalent ring carbons |

**Consequence for the population estimate.** Y1's 2.85% (171/6000) came from a pure-geometry
pre-filter — "N within 2.6 Å of a transition metal with exactly 1 H and 2 C neighbours" — which
does not test symmetry-distinctness. It is therefore an **over-count** of the motif that is
actually stereogenic. See §7 for the measured figure.

## 6. Results

Encoder behaviour, lever ON (all measured; `mirror` is the z-reflection):

| molecule | lever OFF | lever ON |
|---|---|---|
| `POJJOP` (1 locked amine, sole stereocentre) | mirror byte-identical | `[N@@H]{0}` vs `[N@H]{0}`, `mirror == swap(base)` |
| `RIFGUJ_comp_2` (3 Cu-bound amines) | mirror byte-identical | all three invert, nothing else moves |
| `TEGFET`, `KANYUW`, `LARYEI` | already diverged via a second stereocentre | still diverge, now with explicit amine descriptors |
| `JUCCUH` | pendant `[N@@H]` already captured | **byte-identical** — no eligible donor |

Three-property test (`tests/unit/test_locked_donor.py`):

| property | POJJOP | RIFGUJ_comp_2 |
|---|---|---|
| invariant under input renumbering | **whole string**, 3/3 | **descriptors**, 4/4 (see below) |
| invariant under proper rotation | **whole string**, 3/3 | **whole string**, 3/3 |
| flips under reflection | `mirror == swap(base)` | `mirror == swap(only the [N@] tags)` |

Two scoping notes, both forced by measurement rather than convenience:

* `mirror == swap(EVERY tag)` is **false** on RIFGUJ, correctly. Its 1,3,5-trisubstituted
  cyclohexane ring carbons encode a *relative* (all-cis) arrangement that does not invert under
  reflection — verified with the lever **off**, so it is not something this lane introduced.
  Lane 8 recorded the same shape of problem for an η-arene whose mirror differed only in the
  winding character. The assertion used is still a whole-string comparison; tag *counts* are
  never compared, since a symmetric swap defeats that.
* Full-string renumbering stability cannot be asserted on RIFGUJ. With the lever **off** it
  drifts in 3 of 3 renumberings — ring-carbon tags flip, one loses its tag entirely, and the
  `{2}`/`{3}` slot numbers swap. That is the pre-existing 13% stereo-flip class (Lane 8) plus
  slot drift (Lane 2). Asserting it here would import two other lanes' open defects into this
  lane's guard.

## 7. Known residual: internally compensated donor pairs

`ABIFAV_comp_0` (Pd, two adamantanecarboxylate O donors and two *N*-isopropylcyclohexylamine
donors) emits both descriptors — `{1}` as `[N@@H]` and `{3}` as `[N@H]` — and yet base and mirror
remain byte-identical with the lever ON. The two amine ligands are constitutionally identical and
carry **opposite** configurations, so the pair is internally compensated: reflecting the molecule
exchanges which ligand is which, and the descriptor pair is reflection-invariant. On that reading
the identical encoding is correct and the rigid-mirror oracle's "distinct, RMSD 3.30 Å" verdict is
the known conformation-inflated artifact for flexible ligands (the Y1 overview flags exactly this
for the rigid-mirror chirality scan). `LARYEI_comp_0` shows the same pattern at its two
ethylenediamine-derived nitrogens, and there the string still differs because its backbone carbons
do invert.

**What is not established:** whether the fold in `ABIFAV` is purely the (R,S) compensation
described above, or whether the reflection-sensitivity of *slot assignment* is also contributing.
Distinguishing the two needs a conformer-independent isomer oracle, which this lane did not build.
Recorded rather than resolved.

## 8. Scope explicitly not taken

* **The xfail does not flip.** `test_metal_bound_amine_should_diverge_in_raw_string` remains an
  expected failure, and that is now a statement about the *default*: the fix is gated, and
  "levers-OFF output is byte-identical" is the harder acceptance requirement. The capability is
  guarded permanently in `tests/unit/test_locked_donor.py` instead. Flip the marker when the
  lever is promoted.
* **Generator support.** Not attempted. RDKit drops `[N@]` on the parse the generator performs,
  so selection-by-descriptor (the `_axial_narrow` precedent) would need the descriptor carried
  out of band into the parsed structure first. This is the concrete blocker to promoting the
  lever, and it is a generator-side task.
* **Trivalent P in practice.** The mechanism covers metal-locked P by construction and the
  notation demonstrably survives for it, but no corpus case was found where a P donor's tag was
  cleared *and* the four-distinct-classes gate passed. A PR₃ donor's three substituents are very
  often not all distinct, and when they are, the existing Zone-A lone-pair path usually already
  produces a label. So the P half of this lane is **implemented and tested but not
  corpus-demonstrated** — say so rather than claiming it closed something measured.

## 9. Reproducing

```bash
cd ../oin-v045-lane6 && export PYTHONPATH=$PWD/src:$PWD
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python

# the guards
$V -m unittest tests.unit.test_locked_donor

# the defect and the fix, side by side
$V -m tools.injectivity.twin_collision tests/fixtures/POJJOP.xyz
OIN_EMIT_LOCKED_DONOR=1 $V -m tools.injectivity.twin_collision tests/fixtures/POJJOP.xyz
```
