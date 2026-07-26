# Lane 8 — stable stereo under atom renumbering

**What this lane was for:** an **unplanned** lane, opened mid-wave after a probe found that
**13% of molecules emit a different absolute stereochemistry when the same 3D structure is
presented with its atoms listed in a different order** — a soundness defect, not a cosmetic one —
and closed by re-deriving fragment chiral tags from the parent geometry behind
`OIN_STABLE_STEREO`, since promoted to default-ON.

Primary sources: `docs/RENUMBERING_INSTABILITY_v0.4.5.md`,
`src/oinsmiles/oin/stable_stereo.py`, `tests/unit/test_stable_stereo.py`,
`tests/unit/test_stable_stereo_mirror.py`, `docs/PROMOTION_GATE_v0.4.5.md`,
`docs/V045_STATUS_2026-07-25.md`.

---

## ELI5

An `.xyz` file is just a list of atoms with their positions — "carbon here, hydrogen there" —
and **the order of the lines carries no chemical meaning whatsoever.** Shuffle the lines and you
have written down the *identical* molecule sitting in the *identical* place; the picture on the
screen would not move by one pixel. So the 1D string this project produces from that file must
come out byte-for-byte identical. It did not: for about 1 in 8 molecules, shuffling the lines
made the string report the **opposite handedness** — the difference between a left hand and a
right hand, which in chemistry can be the difference between a medicine and a poison. That is not
a typo in the output; it means that on any *one* line ordering the answer was right only by luck.
The cause turned out to be a bookkeeping subtlety: the string does not store "this centre is
left-handed", it stores "this centre is left-handed *relative to the order I happen to be listing
its neighbours in*" — and one step of the code rebuilt each ligand and quietly changed that order.
The fix stops translating the old bookkeeping and instead re-reads the handedness straight off the
coordinates, which no line ordering can change.

---

## The work, visually

```
PART 1 — WHY A CHIRAL TAG IS NOT AN ABSOLUTE STATEMENT

An RDKit chiral tag (CHI_TETRAHEDRAL_CW / _CCW; written @@ / @ in SMILES) is a PARITY
relative to THE ORDER THE NEIGHBOURS APPEAR ON THAT ATOM. It is a permutation sign,
not a configuration.

  same 3D centre, two neighbour orders, two different tags — BOTH CORRECT:

        neighbour order (a, b, c, d)              neighbour order (b, a, c, d)
                 b                                         a
                 |                                         |
        a ——— (C*) ——— c   reads "@"          b ——— (C*) ——— c   reads "@@"
                 |                                         |
                 d                                         d
                                               ^^^^ one transposition ⇒ parity flips

  ⇒ a tag is only meaningful TOGETHER WITH the ordering it was recorded against.


PART 2 — WHERE THE REBUILD BREAKS IT   (get_oin_string, src/oinsmiles/utils/xyz2mol.py)

  get_oin_string does NOT use Chem.PathToSubmol or RWMol.RemoveAtom.
  It REBUILDS every ligand fragment from scratch, in three steps:

   PARENT mol (perceived from the xyz, has a 3D conformer)
        │
        │  step 1:  mw.AddAtom(atom)  for each heavy atom
        │           └── copies the parent's chiral tag  V E R B A T I M
        │               (the tag travels; the ordering it referred to does NOT)
        │
        │  step 2:  hydrogens are DROPPED and folded into SetNumExplicitHs
        │           └── ✗ CHANGES the neighbour list on every centre that had an H
        │               ── and WHICH atoms are hydrogens is a function of input numbering
        │
        │  step 3:  bonds re-added by iterating heavy_indices in ASCENDING PARENT INDEX,
        │           taking each neighbour with nbr_idx > old_idx
        │           └── ✗ CHANGES the neighbour ORDER
        │               ── and "ascending parent index" IS the input numbering
        ▼
   FRAGMENT RWMol  ── tag from step 1 is now a parity against an ordering that
                      NO LONGER EXISTS.  If the two orderings differ in parity:

        ordering A  ──► emits  [C@H]        ┐  same coordinates
        ordering B  ──► emits  [C@@H]       ┘  same molecule       ⇒ SOUNDNESS DEFECT
                                               at least one is WRONG


PART 3 — THE FIX (OIN_STABLE_STEREO, src/oinsmiles/oin/stable_stereo.py)

   DO NOT TRANSLATE THE TAG.  RE-DERIVE IT.

   restamp_fragment_chirality(frag_rw, frag_to_parent, parent_conf)
        │
        ├─ collect `prior` = atoms ALREADY carrying CHI_TETRAHEDRAL_CW/CCW
        │     └── SCOPE GUARD: no other atom is touched, so the lever can never
        │         ADD or REMOVE stereochemistry.  tags_in == tags_out.
        ├─ copy the PARENT conformer's coordinates onto a COPY of the fragment
        ├─ Chem.AssignAtomChiralTagsFromStructure(probe)
        │     └── stamps the tag against the FRAGMENT's OWN neighbour order
        ├─ RDKit DECLINED?  (degree-3, zero-H heteroatom — i.e. EVERY phosphine or
        │  amine donor once the metal bond is cut)
        │     └── _tag_from_geometry(): sign of the triple product of the unit
        │         vectors to the first three neighbours IN BOND ORDER;
        │         positive ⇒ CCW.  |t| < _PLANARITY_TOL (0.05) ⇒ DECLINE
        └─ still nothing ⇒ KEEP the existing tag
              (an absent tag is a worse answer than an uncertain one)

   GEOMETRY IS NOT A FUNCTION OF ATOM NUMBERING  ⇒  order-independent BY CONSTRUCTION


PART 4 — WHY STABILITY ALONE WOULD BE WORTHLESS  (the paired assertions)

        ┌──────────────────────────────────────────────────────────────────┐
        │  A CONSTANT descriptor is PERFECTLY STABLE and encodes NOTHING.  │
        │  The Y2 wave shipped exactly that: sorting symmetry-equivalent   │
        │  axes by sign made an axial token REFLECTION-INVARIANT, and      │
        │  every guard written against the one easy fixture PASSED.        │
        └──────────────────────────────────────────────────────────────────┘

   so each stability assertion is PAIRED with a faithfulness assertion:

     STABLE     renumber → byte-identical            (must hold for EVERY seed)
       ⊕
     FLIPS      z-mirror → DIFFERENT string          (must not collapse)
       ⊕
     TRUE       matches the absolute configuration the fixture is NAMED for
       ⊕
     COUNTED    tag COUNT preserved  ⇒ stability was not bought by emitting fewer tags

LEGEND
  ✗              the step that breaks the invariant
  @ / @@         SMILES chirality tags (CHI_TETRAHEDRAL_CCW / CW)
  parity         sign of the permutation of an atom's neighbour list
  lever          an OIN_* environment variable read through oin/levers.py::lever_enabled
  ⊕              "and, inseparably"
```

---

## Initial assumptions and hypothesis

This lane was **not on the v0.4.5 plan**. v0.4.5's thesis was that the *emitted* string should be
canonical, where previously only the internal comparison key was — a canonicality problem. The
prevailing assumption was therefore that any remaining drift was cosmetic: slot renumbering,
fragment ordering, aromaticity spelling.

Three specific beliefs were held going in, and the third was **written down and handed to the Lane
8 agent as a briefing**:

1. The encoder was assumed **presentation-invariant** in the ways that matter. Nothing had ever
   tested "does the encoder consistently report *one* enantiomer for one structure?"
2. Any renumbering drift was assumed to be **canonicality**, in the same class as the slot and
   body drift the other lanes were closing.
3. **The 13% stereo-flip class was assumed to be downstream of order-dependent bond-order
   perception**, and therefore expected to be fixed for free by Lane 1's
   `OIN_CANONICAL_PERCEPTION`. This was reasoned out and communicated in writing. **It was
   measured and it is WRONG** — see "Dead ends".

The probe's design assumption, which held perfectly, was that **holding the molecular graph fixed
makes the expected answer byte-identical rather than merely similar** — so the ground truth is
known, not inferred.

---

## What was actually found

### The measurement (`tools/canonicality_probe.py`)

Run at `main` @ `20044883`, **all levers OFF**. Sample: **225 molecules (seed 42)** from the
**25,197-basename** tmCAT-tmPHOTO corpus, **2 trials per transform**; **223 encoded
successfully**.

Three transforms, each holding the molecular graph fixed:

| transform | what changes |
|---|---|
| `rotate` | a random **proper** rotation (det = +1) of the coordinates |
| `renumber` | the order atoms appear in the XYZ file |
| `both` | renumber, then rotate |

Rotations are forced **proper** because an improper operation mirrors the structure, which
legitimately changes a chiral molecule's encoding.

### CONFIRMED — headline results

| | count | share |
|---|---:|---:|
| byte-stable across all transforms | 125 / 223 | **56.1%** |
| drifted | 98 / 223 | **43.9%** |
| of which: comparison **KEY** also changed | 47 / 223 | **21.1%** |

Drift by transform: `renumber` **132**, `both` **123**, **`rotate` 0**.

### CONFIRMED — the encoder is fully orientation-invariant; every defect is atom-numbering

Across 225 molecules × 2 trials, **not one** structure changed its OIN string under a random
proper rotation. `_align_to_pai` does its job for orientation. **Every** defect below is an
atom-numbering dependence.

### CONFIRMED — the severity split under pure renumbering, and the 13% is a soundness class

| class | count | share |
|---|---:|---:|
| **stereo-only flip** — string identical after deleting `@`/`@@`, tags differ | **29** | **13.0%** |
| other drift — skeleton, slot numbering, or aromaticity also differs | 79 | 35.4% |
| geometry classification changed (`[M_XXX]`) | 1 | 0.4% |
| aromaticity perception collapsed (> 2 atoms) | 8 | 3.6% |

The stereo-only class is measured **strictly**: the two strings are byte-identical once chiral
tags are removed. **Because SMILES chirality is defined relative to the order neighbours appear
*in the string*, and that order is identical here, the differing tags denote genuinely different
configurations. At least one of the two is therefore wrong.** That is what makes this a soundness
defect and not a formatting one.

Worked example, `FEQFIS_comp_0` — verified the two inputs are a **pure permutation** (identical
sorted coordinate multiset, so no geometry change and no mirroring):

```
BASE : [Au_LIN].C[C@@H](c1ccccc1)N([C@@H](C)c1ccccc1)p{0}1oc2c(...)...[Cl]{1}
RENUM: [Au_LIN].C[C@H](c1ccccc1)N([C@@H](C)c1ccccc1)p{0}1oc2c(...)...[Cl]{1}
                    ^^^^^
```

### CONFIRMED — the root cause is inside the fragment rebuild, not in perception

**Measured: the parent mol's perceived graph and chiral tags were IDENTICAL between the two
orderings for 26 of the 29 affected molecules** — the divergence appeared for the first time
*inside* `get_oin_string`'s rebuild. This is the single measurement that localises the defect, and
it is what refuted the perception hypothesis.

The mechanism, exactly: `RWMol.AddAtom` copies the parent's chiral tag **verbatim**, but a chiral
tag is a **parity relative to the neighbour order on that atom** — and the rebuild changes that
order **twice**, both times as a function of the input atom numbering: hydrogens are folded into
`SetNumExplicitHs`, and bonds are re-added by **ascending parent index**.

### CONFIRMED — the lever closes the stereo class

`OIN_STABLE_STEREO` on `swimlane/v045-lane8`, measured over **10 of the 29** known stereo-flip
molecules, **3 trials × 3 transforms**:

| | byte-stable | drifted | key-level defects |
|---|---|---|---|
| lever OFF | **0/10** | 10/10 | **10** |
| `OIN_STABLE_STEREO=1` | **8/10** | 2/10 | **1** |

### CONFIRMED — and this is the guard that mattered — the descriptor is still faithful

Run precisely **because Lane 8's own WIP never confirmed it**. Tested by z-mirroring each molecule
with the lever ON:

- **10/10 mirrors produce a DIFFERENT string** — nothing collapsed.
- On every case inspected, **`mirror == swap(base)`**, where `swap` exchanges `@` ↔ `@@`: the
  mirror string is the base string with *every* chiral tag inverted and **nothing else changed**.
  That is a textbook enantiomer pair, i.e. the descriptor is doing exactly its job.

### CONFIRMED — the promotion evidence (all six canonicality levers together)

From `docs/PROMOTION_GATE_v0.4.5.md`, measured on `trial/v045-integration` with
`tools/canonicality_probe.py --n 300 --trials 2`, seed 42 fixed so every arm samples the same
molecules; **298 / 299 molecules encoded**:

| arm | byte-stable | comparison **key** broken |
|---|---|---|
| all levers OFF | 173/298 — **58.1%** | 60 — **20.1%** |
| all levers ON | 208/299 — **69.6%** | 16 — **5.4%** |
| **delta** | **+11.5 pts (+35 molecules)** | **60 → 16, a 73% reduction** |

`rotate` drift is **0 in both arms**. What the six promoted levers have in common, and why it made
them safe: **each repairs a renumbered presentation without rewriting the canonical answer.**

### CONFIRMED — the historical-numbers consequence

The comparison key is the harness's acceptance predicate and the basis of every accuracy figure
this project reports. It changed under renumbering for **21.1%** of sampled molecules, so
**reported round-trip accuracy carries a systematic error term that has never been accounted
for.** This does not invalidate relative A/B results (both arms share the input ordering), but
absolute accuracy figures must be read with it in mind. Stated more sharply in the status doc:
**"100% round-trip accuracy" is not currently well-defined**, because round-trip success is partly
a property of the input file's atom order.

### CONFIRMED — why neither existing instrument could see this

- **The Y1/Y2/Y3 injectivity audit** asked *"does the encoder **separate** two enantiomers?"* —
  mirror-twin collision probes. It never asked *"does the encoder **consistently report** one
  enantiomer?"* A mirror-twin probe compares two **different** structures; this probe compares
  **one structure with itself.**
- **The round-trip sweep** compares an input encoding against a *generated* structure's encoding.
  The generator builds coordinates from the OIN string, so **both sides inherit whatever
  configuration the encoder chose on that particular input ordering** — the error is common-mode
  and cancels. Only re-presenting the *same* input differently exposes it.

---

## What was done

### The lever and its wiring

| item | value |
|---|---|
| lever name | **`OIN_STABLE_STEREO`** |
| default | **ON** since v0.4.5 (registered in `src/oinsmiles/oin/levers.py::_DEFAULT_ON`) |
| implementation | `src/oinsmiles/oin/stable_stereo.py` (163 lines) |
| entry point | `restamp_fragment_chirality(frag_rw, frag_to_parent, parent_conf) -> int` (count of tags flipped) |
| helper | `_tag_from_geometry(mol, conf, idx)` |
| constants | `_TETRAHEDRAL = (CHI_TETRAHEDRAL_CW, CHI_TETRAHEDRAL_CCW)`; `_PLANARITY_TOL = 0.05` |
| hook | `src/oinsmiles/utils/xyz2mol.py:1625` — `if lever_enabled("OIN_STABLE_STEREO") and not is_metal and mol.GetNumConformers():` (a 14-line block inside the per-fragment rebuild, placed after the E/Z stereo-atom carry and before `frag_mol = mw.GetMol()`) |

### `restamp_fragment_chirality`, decision by decision

1. **Only atoms that already carry `CHI_TETRAHEDRAL_CW`/`CCW` are restamped**, and no other atom
   is touched. This is the scope guard that makes the lever incapable of adding or removing
   stereochemistry: **the number of tagged atoms going in equals the number going out.** Making
   the string stable by *dropping* a tag is the failure mode this lane exists to prevent.
2. **Bail out fail-safe** if `prior` is empty, or if `frag_to_parent` does not cover every fragment
   atom.
3. **Work on a copy.** `Chem.RWMol(frag_rw)`, conformers removed, parent coordinates copied in via
   `frag_to_parent`, then `UpdatePropertyCache(strict=False)` + `FastFindRings` +
   `AssignAtomChiralTagsFromStructure`. A perception failure must never corrupt the fragment that
   is about to be serialised — the whole block is wrapped in `except Exception: return 0`.
4. **RDKit declines a whole class this encoder produces**, so a geometric fallback is required.
   `AssignAtomChiralTagsFromStructure` returns `CHI_UNSPECIFIED` for a **degree-3 atom with no
   hydrogen** — and *after the metal bond is cut, every phosphine or amine donor looks exactly like
   that*. Without the fallback the order-dependent copied tag would survive on precisely the donor
   atoms this notation cares most about. `_tag_from_geometry` therefore reads the handedness
   straight off the coordinates using **RDKit's own convention** (`assignChiralTypesFrom3D`): take
   the first three neighbours **in the atom's own bond order**, sign the triple product of the unit
   vectors to them; the implicit fourth substituent (a hydrogen or a lone pair) sits opposite, and
   a **positive** triple product is `CHI_TETRAHEDRAL_CCW`.
5. **`_PLANARITY_TOL = 0.05`** — below this `|triple product|` of the three *unit* neighbour
   vectors the centre is too close to planar to read a handedness from, so the existing tag is kept
   rather than resolved by numerical noise. **Validated against RDKit on every tagged atom of the
   project's 48 XYZ fixtures: 255 agreements, 3 disagreements, all three at
   `|unit triple product| < 0.016`** — planar centres neither method should be reading, and all
   three were Cp ring carbons at `|t| = 0.011–0.016`.
6. **If nothing can be derived, keep the tag we have.** An absent tag is a worse answer than an
   uncertain one.

**Rejected alternatives:**

- **Translate the parity instead of re-deriving it** — i.e. compute the permutation between the
  parent's and the fragment's neighbour orders and flip the tag when the permutation is odd.
  Rejected: it is a re-implementation of the exact bookkeeping that broke, and it would have to be
  right for three separate order-changing steps (H folding, bond re-add order, and any later
  sanitisation). Geometry is not a function of atom numbering, so re-deriving is
  order-independent **by construction** rather than by care.
- **Extract fragments with `Chem.PathToSubmol` / `RWMol.RemoveAtom` instead of rebuilding.**
  Not taken: the rebuild is load-bearing for a great deal else in `get_oin_string` (hydrogen
  folding, slot markers, E/Z carry), and replacing it is a far larger blast radius than
  re-stamping four to eight atoms per fragment.
- **Narrow `core/chirality.py`'s Zone-A `total_degree < 4` clear** so metal-bound donors keep
  their tags. Explicitly out of scope: the restamp only *corrects tags that already exist*, so
  making metal-locked donor tags **exist in the first place** is Lane 6's job (P3). The status doc
  records the ownership call as *"Do not widen Lane 8."*

### Guard tests — the paired-assertion design

`tests/unit/test_stable_stereo.py`

Helpers: `_encode(path, stable)` sets the lever **explicitly in both directions** (see the trap
below); `_renumbered(path, tmpdir, seed)` shuffles atom lines carrying coordinates across
untouched; `_mirrored(path, tmpdir)` reflects through the yz plane; `_cip_labels(oin)` strips
`{slot}` markers and the `_GEO` suffix, parses each fragment with `MolFromSmiles`, and labels with
**`rdCIPLabeler`** (the rigorous implementation, not legacy `AssignStereochemistry`).

Fixtures:

- `_UNSTABLE = "ROGYAO_comp_0.xyz"` — 49 atoms, one Rh, tetrahedral carbon stereocentres on a
  bidentate ligand. Chosen from the 29 renumbering-unstable molecules because it is small,
  encodes quickly, and **drifts in the stereo tags WITHOUT any accompanying skeleton or slot
  drift** — so a failure there can only be this defect.
- `_RR_FIXTURES = ("PdCl2-RR-BDPP.xyz", "PdCl2-RR-BDNN.xyz")` — (2R,4R)-pentane-2,4-diyl
  backbones. **Their absolute configuration is checkable by hand from the raw coordinates**:
  priority at each centre is D > CH2 > CH3 > H with D = P or N, and the triple product of the unit
  vectors to priorities 1,2,3 is negative for both centres in both fixtures, so **the correct
  answer is (R,R), independently of anything RDKit or this codebase computes.**

| class / test | what it asserts |
|---|---|
| `TestStereoIsStableUnderRenumbering::test_lever_off_reproduces_the_defect` | **Anti-vacuity guard.** The fixture must really be unstable, else the stability assertion below proves nothing. Asserts over the **whole** seed set, not any particular seed |
| `::test_lever_on_is_byte_stable_under_renumbering` | byte-identical for **every** seed in `SEEDS` |
| `::test_lever_on_keeps_reviewed_complexes_stable` | BDPP, BDNN, `Rh-RR-DIPAMP-Cl2` at seeds 4 and 9 |
| `TestStableStereoIsCorrect::test_lever_on_emits_the_true_configuration` | the (2R,4R) fixtures encode as `["R","R"]` |
| `::test_lever_off_emits_the_inverted_configuration` | documents the defect: the unpatched path emits `["S","S"]` |
| `::test_mirror_image_inverts_every_tag` | constitution survives, `base != mirror`, tag counts equal, **and every individual tag differs** |
| `::test_achiral_fixture_mirror_is_IDENTICAL` | the **over-sensitivity** guard, in the opposite direction |
| `::test_mirror_of_an_rr_fixture_is_ss` | CIP of the mirrored input is `["S","S"]` |
| `::test_tag_count_is_preserved` | **stability must not be bought by emitting fewer tags** |
| `TestLeverIsOnByDefault::test_unset_env_matches_explicit_ON` | unset takes the promoted default |
| `::test_registry_agrees_that_the_lever_ships_on` | `OIN_STABLE_STEREO in levers.default_on()` |
| `::test_lever_does_not_change_the_ORIGINAL_ordering` | ON == OFF on the input's own ordering — *the desirable property*, see below |

**The seed set is `SEEDS = tuple(range(12))`, deliberately wide.** The two assertions are
**asymmetric**: stability must hold for **EVERY** seed, while "the fixture really is unstable"
needs only **ONE** seed to drift. So widening strengthens the first and de-flakes the second at
once. It was originally `(1, 2, 3, 5, 8)` and that was **too narrow**: none of those five happens
to trigger `ROGYAO_comp_0`'s drift, so `test_lever_off_reproduces_the_defect` failed even though
the molecule genuinely **IS** unstable — `canonicality_probe.py --only ROGYAO_comp_0 --trials 5`,
which draws its own permutations, sees **5 of 15 transform-trials drift** (subclass
`rdkit_canonical`, key-level too). A five-sample miss on a ~1-in-3 event is unremarkable;
**hardcoding five seeds was the bug, not the fixture.**

`tests/unit/test_stable_stereo_mirror.py` — the **permanent** mirror guard (commit `92c11616`),
over `EJUJUP_comp_0.xyz` and `OCUGIC_comp_0.xyz`, both from the measured set of 29:

| test | what it asserts |
|---|---|
| `TestStableStereoIsBothStableAndFaithful::test_renumbering_is_byte_stable` | half one, 3 seeds per fixture |
| `::test_mirror_does_not_collapse` | half two, **the one with teeth**: the mirror must differ, **and** once every stereo-bearing token is neutralised the two strings must be identical — so the difference is confined to configuration and did not smuggle in a constitution change |
| `::test_pure_sp3_case_inverts_every_chiral_tag` | the strict form for `OCUGIC_comp_0`: `_swap_tags(base) == mirror` |

---

## Dead ends and refutations

### ⚠ THE BIG ONE — the fix EXPOSED A FOUR-MONTH-OLD WRONG ANSWER, blessed by a CIRCULAR oracle

`tests/unit/test_chiral_p.py` and `tests/unit/test_chiral_n.py` asserted **`["S", "S"]`** for
fixtures named **`PdCl2-RR-BDPP`** and **`PdCl2-RR-BDNN`**, citing *"verified by RDKit CIP"*.

**The verification was circular.** It took the encoder's **OWN emitted string**, reparsed it with
`MolFromSmiles`, and ran CIP on the result. **`rdCIPLabeler` converts a parity tag into an R/S
label — it does not check that tag against anything.** Hand it an inverted tag and it returns an
inverted label with full confidence. So the "oracle" was a snapshot of the encoder's output, and
an inverted tag was self-consistent and passed. (The old form also used the **legacy**
`AssignStereochemistry` rather than `rdCIPLabeler`.)

**Ground truth is the geometry, the one thing no encoder bug can rewrite.**
`Chem.AssignStereochemistryFrom3D` on the **parent complex** gives **R at both centres**, agreeing
with the **`(2R,4R)`** in the fixtures' own filenames. Therefore **`OIN_STABLE_STEREO` is CORRECT
and the goldens recorded the DEFECT.**

Repaired in commit **`8f715699`**. Both test modules now derive truth from coordinates **AND**
cross-check that the emitted string agrees — **two tests, not one**, which is what closes the loop
the circular form left open:

- `test_chiral_p.py::TestChiralP::test_p_cip_from_geometry` — `get_tmc_mol(..., with_stereo=True)`
  → `AssignStereochemistryFrom3D` → `rdCIPLabeler` → `["R","R"]`, with the failure message
  spelling out *"if this fails the fixture or the perception changed, NOT the golden string — do
  not 'fix' it by editing the expectation."*
- `test_chiral_p.py::TestChiralP::test_emitted_string_agrees_with_the_geometry` — reads the emitted
  ligand body back and requires it to agree with the independent geometric answer.
- The same two, `test_n_cip_from_geometry` / `test_emitted_string_agrees_with_the_geometry`, in
  `test_chiral_n.py`.
- The golden was re-pinned (both backbone tags flipped `@@`/`@` → `@`/`@@`), with the reason in the
  module docstring: **both spellings describe R at both centres**; the per-atom symbol differs from
  the label because parity is relative to neighbour write order.

New golden for BDPP:

```
[Pd_SPL].C[C@H](C[C@@H](C)P{0}(c1ccccc1)c1ccccc1)P{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}
```

Note the corroborating detail from the triage: **BDPP/BDNN round-trip keys DO change, correctly**
— unlike CisPlatin, TransPlatin, `fac-Ir(ppy)3` and BINAP, which relabel under
`OIN_CANONICAL_SLOTS` with `canonical_roundtrip_key` **identical** (relabelings, not new isomers).
The BDPP/BDNN key change is the signature of the old tags having been inverted.

### ⚠ A vacuous guard that was GREEN — "unset means off"

`test_stable_stereo.py::TestLeverIsOffByDefault` asserted **"unset env == explicit OFF"**. That
premise stopped being true the moment the lever was promoted to default-ON — **yet the test kept
PASSING**, because on the fixture's **original atom ordering the lever is a no-op** (which is
exactly what `test_lever_does_not_change_the_ORIGINAL_ordering` now pins). **A guard asserting a
false premise, green.**

Renamed and inverted to `TestLeverIsOnByDefault`, asserting the lever ships **ON**. It also gained
`test_registry_agrees_that_the_lever_ships_on`, so the claim is checked against
`levers.default_on()` and not just against behaviour.

Found by `tests/unit/test_levers.py::TestNoTestUnsetsAPromotedLever`, which exists because **this
"unset means off" pattern has cost 23 test failures across two promotions** — and this instance was
invisible to all of them. In the v0.4.5 promotion triage the pattern accounted for **17 of the 36
broken tests**: every lever-OFF test expressed "off" by *deleting* the env var, which was correct
while all levers defaulted off and silently means **ON** now. Lane 8's
`test_lever_off_reproduces_the_defect` was, for a while, asserting the defect against the **fixed**
code — it failed loudly, which is the only reason this was caught at all.

### REFUTED — "the stereo class shares Lane 1's bond-order-perception root"

This hypothesis was reasoned out and **told to the Lane 8 agent in writing**. Measured on the three
worked examples, **3 trials × 3 transforms**:

| levers | byte-stable | key-level defects remaining |
|---|---|---|
| all OFF | 0/3 | **3** |
| `OIN_STABLE_METAL_AC` | 0/3 | **2** |
| `+ OIN_CANONICAL_PERCEPTION` | 0/3 | 2 |
| `+ OIN_CANONICAL_BODY + OIN_CANONICAL_SLOTS` (all four) | 0/3 | **2** |

**Only the metal-AC fix closes anything. `FEQFIS_comp_0`'s stereo flip and `CEBVIR_comp_0`'s
aromaticity collapse survive all four levers.** So the stereo class does **not** share a root with
Lane 1's perception work, and Lane 8's independent investigation was genuinely required rather than
redundant. Scope caveat recorded with the table: **n = 3, deliberately the hardest hand-picked
cases.** It does not refute Lane 1's general result (6 fixed / 0 regressed over 250 molecules); it
refutes only the specific claim that Lane 1's lever closes the stereo-flip class.

Practical consequence: the lane was told to **finish** its WIP (`swimlane/v045-lane8` @
`8fdccb55`), not fold it into Lane 1. Note the lane had already merged Lane 1's lever
(`220c191b`) purely so the shared-root question could be *tested*; that merge is behaviourally
inert with no env set.

### REFUTED — "counting tags checks that the mirror still differs"

An earlier hand-run of the mirror check reported **7/10** by comparing tag **counts**. **That test
is blind to a symmetric swap** — a molecule with three `@@` and three `@` maps to three `@` and
three `@@`, identical counts. The 7/10 was an **artifact of the test, not a property of the fix**;
comparing the whole string under an `@` ↔ `@@` swap gives **10/10**. Both permanent guards now
compare whole strings.

### REFUTED — "mirror == swap(@ tags)" as a universal assertion

**This notation carries stereochemistry in TWO places.** `EJUJUP_comp_0` is a Cr arene whose mirror
differs from the original **solely in the eta winding character, `{0>}` → `{0<}`** — the
coordinated ring face — with its `[C@@H]` tag **unchanged**. So the strict swap form is wrong in
general and fails on eta complexes **for a correct fix**. `test_mirror_does_not_collapse` therefore
asserts *where* the difference lives (neutralise `@`, `@@`, `>` and `<`, then require equality)
rather than asserting a specific tag flip; the strict form is kept only for a pure-sp3 fixture and
self-skips if it sees eta winding.

### REFUTED — `ROGYAO_comp_0` is chiral

It carries two `[C@@H]` tags and reads like an obvious enantiomer test, so
`test_mirror_image_inverts_every_tag` originally demanded its mirror differ. **It does not, and it
should not: Lane 7's torsion-aware oracle returns `RIGID_ACHIRAL` — the mirror superimposes at
0.423 Å over 4 automorphisms, against a 0.5 Å threshold.** Removed from that test and replaced by
`test_achiral_fixture_mirror_is_IDENTICAL`, asserting the **opposite** — the over-sensitivity guard
that a naive "mirrors must always differ" rule gets backwards. An encoder that distinguished this
molecule would be **manufacturing stereochemistry that is not there**, a false positive the round
trip would then fail on forever.

> ⚠ **Why it looked chiral, and the instrument lesson.** The **older rigid** oracle
> (`tools/injectivity/oracle.py`) reported *"distinct, mirror RMSD 2.586 Å, ENCODER-BLIND
> (total)"* — a **false positive from hitting its 4000-automorphism cap**. It enumerates on the
> **H-explicit** graph, where this molecule's four methyls starve the budget, inflating rigid RMSD
> on methyl-rich species. Lane 7 found and documented that instrument defect and deliberately left
> it in place; **this is the concrete case where it would have sent someone chasing a bug that does
> not exist. Use `torsion_oracle.py`.**

### Not a dead end, but a WIP claim that later measurement OVERTURNED

The Lane 8 WIP commit (`8fdccb55`) recorded, in the agent's own last words, that *"the fix leaves
lone-pair P centres unhandled — RDKit declines degree-3/0-H"*, and named trivalent phosphorus
donors as its residual. Two later findings reframe that:

1. **Ownership.** `stable_stereo.py:112-118` only corrects tags that **already exist**, so a P
   donor whose tag was cleared by the Zone-A rule (`core/chirality.py:722-727`) has nothing to
   restamp. Making metal-locked donor tags exist is **Lane 6 (P3)**. The two lanes are
   complementary, not overlapping.
2. **The gap itself is smaller than briefed.** Lane 6 measured over **400 corpus molecules: all 7
   eligible metal-locked P donors already carried a tag** from the existing Zone-A lone-pair path
   — the restore was a no-op for every one. Worse for the original framing, **the two molecules
   cited as Lane 8's residuals are not stereocentres at all**: `FEQFIS`'s metal-bound P is aromatic
   P(N)(O)(O)Au with only 3 distinct symmetry classes (the two oxygens are equivalent), and
   `CEBVIR`'s N donors are pyridine-type with 3 neighbours. **Emitting for either would be the
   over-sensitivity failure, not a fix.** The status doc's instruction is explicit: *"Do not cite
   this lane as having closed a phosphorus defect."*

**Still open regardless:** 2 of 10 molecules keep drifting, **1 at key level.**

---

## Where it landed

**Branch `swimlane/v045-lane8`, tip `b7355cfa`. Fully merged** — `main` is 162 commits ahead and 0
behind, so `git log main..swimlane/v045-lane8` is empty. Merged via `05e688cf`
(→ `trial/v045-integration`) and `6eb82071` (→ `release/v0.4.5`).

Commits, in order:

| commit | subject |
|---|---|
| `220c191b` | `merge(lane8): take Lane 1 canonical-perception lever for shared-root testing` (behaviourally inert with no env set) |
| `8fdccb55` | `WIP(lane8): OIN_STABLE_STEREO scaffold, INCOMPLETE — agent halted mid-task` (preserved by the orchestrator after all lane agents were terminated by a spend limit; **not** a finished lane) |
| `92c11616` | `test(lane8): permanent mirror guard for OIN_STABLE_STEREO` |
| `b7355cfa` | `fix(lane8): repair test_stable_stereo -- missing fixture and three wrong assertions` |

Landed **outside** the lane branch but load-bearing for this lane:

| commit | what |
|---|---|
| `1450b5ce` | `release(v0.4.5): integrate 16 lanes and PROMOTE the six canonicality levers` — **Wave D**; `OIN_STABLE_STEREO` promoted to default-ON via the new single-source registry `src/oinsmiles/oin/levers.py` |
| `8f715699` | `fix(v0.4.5): triage the 36 promotion failures — 2 real defects, 1 wrong golden` — the circular-CIP-golden repair, the "unset means off" repair across 17 tests, and the `TestLeverIsOffByDefault` → `TestLeverIsOnByDefault` inversion |

Final state:

| item | value |
|---|---|
| lever | `OIN_STABLE_STEREO` |
| default | **ON** (`levers.py::_DEFAULT_ON`, alongside `OIN_BORON_CAGE`, `OIN_CANONICAL_BODY`, `OIN_CANONICAL_PERCEPTION`, `OIN_CANONICAL_SLOTS`, `OIN_CANONICAL_ETA_WINDING`, `OIN_STABLE_METAL_AC`) |
| opt out | `OIN_STABLE_STEREO=0` — and it now genuinely means off; `_FALSEY = {"0", "", "false", "no", "off"}` |
| src files | `src/oinsmiles/oin/stable_stereo.py` (new), `src/oinsmiles/utils/xyz2mol.py` (+14-line hook at `:1625`) |
| guard modules | `tests/unit/test_stable_stereo.py`, `tests/unit/test_stable_stereo_mirror.py` |
| fixtures added | `tests/fixtures/EJUJUP_comp_0.xyz`, `tests/fixtures/OCUGIC_comp_0.xyz` (`8fdccb55`), `tests/fixtures/ROGYAO_comp_0.xyz` (`b7355cfa`) |
| suite | `b7355cfa` closed **5 errors** in the integrated suite (729 tests, 717 OK, 5 errors, **all one root cause**: `ROGYAO_comp_0.xyz` was referenced but never `git add`ed — it existed only in the gitignored dataset, so every test in the file errored with `FileNotFoundError`). These two modules are now **14 OK**. |

Reproduce:

```bash
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py --n 300 --trials 2 \
    --out tmCAT-tmPHOTO_xyz_dataset/results-v0.4.5-canonicality/baseline-main
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py \
    --only FEQFIS_comp_0,DUDREA_comp_0,CEBVIR_comp_0 --trials 3 -v
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py --only ROGYAO_comp_0 --trials 5
```

> ⚠ **Stale docstrings to fix on a clean tree.** `tests/unit/test_stable_stereo.py`'s module
> docstring still says *"`OIN_STABLE_STEREO` (default OFF)"*, and the in-code comment at
> `xyz2mol.py:1625` still says *"Default OFF -> byte-identical output."* Both predate the Wave-D
> promotion and are now wrong. `TestLeverIsOnByDefault`'s own docstring dates the rename to
> v0.4.6 while the promotion it responds to landed in the v0.4.5 release commit — read the commit,
> not the docstring, for the date.

---

## Open questions / for the next agent

1. **The 2 residual molecules of 10, one at key level, are not explained.** The trivalent-P story
   that was offered as the explanation has since been measured and does not hold (all 7 eligible
   metal-locked P donors already carried tags; the two cited molecules are not stereocentres).
   The residual therefore has **no current diagnosis**.
2. **The 8 aromaticity-collapse molecules (3.6%) are untouched.** `CEBVIR_comp_0` survives all
   four canonicality levers. Suspect named in the doc: `AC2BO`'s arbitrary resonance form and
   `get_UA_pairs`' non-unique `nx.max_weight_matching` being atom-order dependent
   (`utils/xyz2mol_local.py:800`, `:542`).
3. **Only 10 of the 29 known stereo-flip molecules were measured.** A full 29-molecule arm was
   never run; the promotion decision rests on the 300-molecule aggregate instead.
4. **Absolute accuracy figures across the project's history still carry an unaccounted error
   term** (21.1% key instability at levers-OFF, 5.4% with all six on). Any restatement of a
   historical pass rate should say which lever configuration produced it.
5. **The circular-oracle pattern is a class, not an incident.** The rule to carry forward:
   **`rdCIPLabeler` RELABELS a tag; it never CHECKS one.** Any assertion about absolute
   configuration must be anchored on `AssignStereochemistryFrom3D` over the parent geometry (or on
   hand-computed triple products), and the string check must be a *second*, separate test. Audit
   any remaining test that runs CIP on reparsed encoder output.
6. **The "unset means off" trap now has a lint** (`test_levers.py::TestNoTestUnsetsAPromotedLever`)
   and an in-test marker convention (`# lever-lint: intentional-unset`). It has cost 23 failures
   across two promotions. Do not add a lever-OFF test that expresses "off" by deleting the
   variable.
7. **`OIN_EMIT_LOCKED_DONOR` (P3) is built, validated, and unusable in the shipped default**,
   because `OIN_CANONICAL_BODY` (default-ON) reparses the body and sanitising the metal-free
   fragment clears `[N@]` on a 2° amine. `levers.py::_HELD_OFF` records that the obvious fix was
   **tried in v0.4.6 and MEASURED wrong** — setting the tag after the sanitise moves the canonical
   write order, and `@`/`@@` is a parity relative to that order. This is the nearest live neighbour
   of Lane 8's mechanism and the same lesson in a different place.
