# Lane 5 — Metal-centred Δ/Λ helicity (Y1 blind spot **P1**)

**The blind spot:** a tris(bidentate) metal complex comes in two mirror-image forms, Δ and Λ,
and the encoder wrote no descriptor at all for the metal centre — so the two enantiomers
collapsed onto one OIN at the round-trip key, silently.

> **This lane was NEVER STARTED in v0.4.5. It was built in v0.4.6.** The v0.4.5 status doc
> records the reason: **0 of 150** corpus molecules emit a metal `@` tag today, so the lane would
> be *creating* a descriptor rather than un-folding a collapse — a larger piece of work than the
> lane budget assumed — and it needed the Δ/Λ tris-bidentate fixture that Lane 7 built.
> `swimlane/v045-lane5` does not exist. The nine Lane 5 commits are on `main`, dated 2026-07-26,
> shipped in release `d799de1f` (v0.4.6).

---

## ELI5

Wrap three identical two-armed ligands around a metal atom and they have to spiral, like the
threads of a screw. A right-handed spiral and a left-handed spiral are **mirror images**: hold
one up to a mirror and you see the other, and no amount of rotating the real one will ever make
it match its reflection — the same way your left hand never fits into a right-handed glove.
Chemists call these Δ and Λ, and they are genuinely different substances. OIN wrote *nothing at
all* about the metal's handedness, so both spirals got the same string; worse, the tool that
checks whether a round trip succeeded deliberately ignored metal handedness, so nobody could
ever have noticed. Building a descriptor for this took **four attempts, three of which
measurement destroyed** — the final answer being that the handedness does not live in where the
donor atoms sit, but in *which donors belong to the same ligand*.

---

## The work, visually

```
 THE FOUR FORMULATIONS, AND WHAT KILLED EACH
 ===========================================

 (1) SIGNED VOLUME of the first four donors in CANONICAL SLOT ORDER
     ┌──────────────────────────────────────────────────────────────────────┐
     │ measured OK: ZUMNEC normalized triple product -0.862                 │
     │              JEGKOW                          +0.026   (30x margin)   │
     │              inverts under reflection; invariant under 6 rotations   │
     └──────────────────────────────────────────────────────────────────────┘
     ✗ KILLED BY: renumbering. ZUMNEC is HOMOLEPTIC -- its six O donors are
       symmetry-equivalent, so element / distance-to-metal / distance-multiset
       ALL TIE. Tie falls to input order. Some resolutions differ by an ODD
       permutation, which INVERTS a signed volume.
             sign flipped  1 -> -1  on 2 of 4 shuffles
       (sub-refutation along the way: an ABSOLUTE volume threshold was wrong in
        KIND -- signed volume scales as (bond length)^3, so no single number
        serves a 2.0 A Rh-N sphere and a 2.7 A Rh-I one. Replaced by the
        dimensionless triple product in [-1,1], band 0.15.)

               │  "stop needing an ordering at all"
               ▼

 (2) MAGNITUDE THRESHOLD on a permutation-invariant OSIPOV-PICKUP-DUNMUR index
     ┌──────────────────────────────────────────────────────────────────────┐
     │ ALL FOUR properties measured and REAL:                               │
     │   6 random donor permutations  -> identical  (was 1 -> -1 before)    │
     │   3 proper rotations           -> identical                          │
     │   reflection                   -> EXACT negation (-4.807e-4/+4.807e-4)│
     │   perfect square (synthetic)   -> exactly +0.000e+00                 │
     │   ideal octahedron (synthetic) -> exactly +0.000e+00                 │
     └──────────────────────────────────────────────────────────────────────┘
     ✗ KILLED BY: real structures. Donors taken from PERCEPTION, not a cutoff:
             ZUMNEC  chiral Delta/Lambda        -4.807e-04
             JEGKOW  ACHIRAL sq. planar, puckered  -3.287e-04
             ideal square (synthetic)           +0.000e+00
       Same order of magnitude, 1.5x apart. Crystallographic PUCKER in an
       achiral complex reads like genuine helicity. NO threshold separates them;
       exact-zero cancellation holds only for IDEALIZED coordinates.
       (and: _CHIRALITY_EPS was first set to 1e-3 -- LARGER than the real
        4.81e-4 signal -- so it would have called ZUMNEC ACHIRAL.)

               │  "chirality is a SYMMETRY property; use a symmetry test"
               ▼

 (3) UNCONSTRAINED MIRROR-SUPERPOSITION symmetry test
     mirror the donor set; ask whether ANY permutation + PROPER rotation superimposes it
     ┌──────────────────────────────────────────────────────────────────────┐
     │ fixture   test says    truth                                        │
     │ ZUMNEC    ACHIRAL      chiral    <-- WRONG                          │
     │ JEGKOW    achiral      achiral                                      │
     │ ideal sq  achiral      achiral                                      │
     └──────────────────────────────────────────────────────────────────────┘
     ✗ KILLED BY: it is RIGHT about the point set and the point set is the
       WRONG INPUT. Six oxygens at octahedral vertices admit improper
       operations -- as a bare point set ZUMNEC IS achiral. A permutation
       search over UNLABELLED points is free to RE-PAIR donors into different
       ligands, which no physical operation can do.
       Retro-explains (2): the -4.807e-04 was residual crystallographic
       distortion. It was NEVER detecting helicity.

               │  "colour the points by which LIGAND they belong to"
               ▼

 (4) CHELATE-AWARE mirror superposition   <<< WHAT WORKS >>>
     permutations restricted to those PRESERVING CHELATE MEMBERSHIP:
     a relabelling may map whole chelates onto whole chelates of the SAME SIZE, nothing else
     ┌──────────────────────────────────────────────────────────────────────┐
     │ ZUMNEC  48 admissible perms   best mirror RMSD 1.3752 A -> CHIRAL   │
     │ JEGKOW   4 admissible perms   best mirror RMSD 0.0582 A -> ACHIRAL  │
     │                                       24x separation                │
     │ _ACHIRAL_RMSD_TOL = 0.35  -- a WIDE BAND, not a tuned constant      │
     └──────────────────────────────────────────────────────────────────────┘
     admissible-permutation counts:  tris-bidentate 3! x (2!)^3 = 48
                                     sizes 1-2-1              =  4
                                     four monodentates  4!    = 24

 THE LESSON
 ==========
   Delta/Lambda helicity is a property of CHELATE CONNECTIVITY, not of donor
   POSITIONS. Reflecting a Delta complex gives Lambda only because the
   reflection cannot be undone WHILE THE CHELATE PAIRING STAYS INTACT.

 THE SHIPPED PIPELINE (all of it default-OFF)
 ============================================
   input .xyz
      │
      ▼  get_tmc_mol()   <- donors = metal's PERCEIVED neighbours (NOT a distance cutoff)
      │                     chelate partition = connected components after deleting the metal
      ▼  token_for_mol(tmc_mol)                          xyz2mol.py:1390-1396
      │    computed BEFORE _align_to_pai, because PAI alignment may REFLECT
      │    the coordinates and a reflection INVERTS a chirality descriptor
      ▼
   [Mo_OCT]....O{5}  |mc:-|            <- trailing sidecar, single return site :2038
      │
      ├──► compare.py::_METAL_CONFIG_TOKEN_RE  ── FOLDS the token out of the key
      │      ⚠ the fold MUST be removed in the commit that promotes the lever
      │
      └──► metallogen_adapter.py: target_mc = parse_metal_config_token(original_oin)
             accept-first requires:  token_for_mol(cmol) == target_mc
             ^ needed SEPARATELY, because the key FOLDS the token, so a key
               match says NOTHING about helicity -- accepting on the key alone
               would hand back the WRONG ENANTIOMER while reporting success

 LEGEND
 ======
   ✗ refuted        ▼ next step in time        ⚠ obligation / trap
   |mc:±|  the opt-in metal-configuration sidecar token
   ZUMNEC  tris(catecholato)Mo -- HOMOLEPTIC, helicity is its SOLE stereogenic element
   JEGKOW  Rh(I) square planar, four DIFFERENT donors (C, P, N, I) -- achiral control
   RMSD    root-mean-square deviation after optimal PROPER-rotation superposition
```

---

## Initial assumptions and hypothesis

### Why the lane was skipped in v0.4.5

`docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` §4 and `docs/agentic-notes/v0.4.5/V045_STATUS_2026-07-25.md`:

> **Lane 5 (metal Δ/Λ, P1) was not started.** 0/150 molecules emit a metal `@` tag today, so it
> would be *creating* a descriptor rather than fixing a collapse — a larger piece of work than
> the lane budget assumed, and it needs the Δ/Λ tris-bidentate fixture Lane 7 built.

That measurement came from Lane 2 (`spec/handoffs/v0.4.5/PROGRESS-lane2.md`): "0 of 150 corpus
molecules emit a metal `@SP1`/`@OH1`". It also retired a Lane 2 risk — metal-`@`-tag staleness
under relabeling is "not a live risk" for the same reason.

### The fixtures Lane 7 built for it (commit `cc0dd117`)

Selected from a full **26,230-structure census**:

- **`tests/fixtures/ZUMNEC.xyz`** — tris(catecholato)Mo. The **only** homoleptic tris-bidentate
  in the corpus whose two O donors per chelate are symmetry-equivalent
  (`CanonicalRankAtoms(breakTies=False)` gives them one rank), so it has **no fac/mer distinction
  at all** and metal helicity is its **sole** stereogenic element. Mirror distinct at
  **1.34 Å**; `@OH` permutation **11 → 9** under reflection.
- **`tests/fixtures/JEGKOW.xyz`** — Rh(I) square planar with four **different** donor elements
  (C-carbonyl, P, N, I). Exercises RDKit's `@SP` path. Its mirror is correctly **not** distinct
  (**0.099 Å**).

**Why this pair closes a trap:** `fac-Ir(ppy)3` — the Y1 fixture — has *unsymmetric* (C,N)
chelates, so a Lane-5 descriptor secretly encoding **fac/mer** rather than **helicity** would
still appear to work on it. ZUMNEC fails such a descriptor. And JEGKOW's non-distinct mirror is
chemically correct rather than a weak fixture: a square-planar complex is planar, so its
coordination plane is a **mirror plane**, and four different donors give **diastereomers, not
enantiomers**. Demanding oracle-distinctness there is a category error; the right distinctness
operator for `@SP` is a **donor swap**, which rides on Lane 7's `swap_donor` twin operator.

### The v0.4.5 handoff's assumption, which turned out wrong in two ways

`spec/handoffs/v0.4.5/Lane5-metal-stereo.md` declared a **hard dependency on Lane 2**: *"it
builds the canonical donor ordering this descriptor is expressed against… Once donors have
canonical slot indices, the metal configuration can be expressed as a permutation of those
canonical slots, reproducible by construction."*

Both halves failed:

1. **No canonical ordering was needed at all** (see formulation 2 — a permutation-invariant sum
   removes the requirement).
2. **Lane 2's ordering would have been the wrong quotient anyway.** Lex-min slot labelling is
   canonical up to the **full automorphism group**, which for a homoleptic sphere includes
   **improper** elements — precisely the elements a chirality descriptor must not quotient by.

The Y2 wave (`docs/agentic-notes/injectivity/INJECTIVITY_Y2_FEASIBILITY.md`) had recorded that RDKit's
`AssignStereochemistryFrom3D` *does* perceive the configuration: `fac-Ir(ppy)3` gets
`CHI_OCTAHEDRAL` with `_chiralPermutation` **10** vs **8** for its mirror. That is a
distinguisher, not a canonical token — RDKit lists non-tetrahedral CIP and canonicalization as
"totally missing".

---

## What was actually found

### Confirmed

| claim | measurement | source |
|---|---|---|
| P1 is **key-blind**, the dangerous cell | `fac-Ir(ppy)3` and its mirror give different **raw** strings, but only by slot renumbering, which `_polyhedron_signature` deliberately folds ⇒ `raw_equal = False`, `key_equal = True` | `docs/agentic-notes/injectivity/INJECTIVITY_Y1_P1_METAL.md` |
| the enantiomers really are distinct | `fac-Ir(ppy)3` min proper-rotation mirror RMSD **3.19 Å** (achiral controls ~0.05 Å); ZUMNEC **1.34 Å** | Y1 doc / `V045_STATUS` |
| the encoder emits no metal stereo at all | `generate_inline_string` builds the metal token as bare `[Metal_GEO]`; the `@…` group in `METAL_REGEX` is a **vestigial reader with no producer** | `oin/inline.py:110`, `:51` |
| the key has *always* folded metal stereo | `_METAL_STEREO_RE` strips `@SP`/`@OH`/`@TB`, `_polyhedron_signature` folds slot relabelling | `oin/compare.py` |
| nothing in the corpus emits a metal tag | **0 / 150** | Lane 2 |
| the OPD index is genuinely permutation-invariant, rotation-invariant, reflection-inverting | 6 random donor permutations identical; 3 proper rotations identical; reflection gives **exact** negation (−4.807e-04 → +4.807e-04) | commit `c72fbc1b` |
| the chelate-aware test separates the fixtures | ZUMNEC **1.3752 Å** (48 perms) vs JEGKOW **0.0582 Å** (4 perms) — **24×** | commit `848af5c6` |
| the shipped emit is correct on all three controls | **verified in this retrospective**: `ZUMNEC → \|mc:-\|`, `ZUMNEC z-mirror → \|mc:+\|`, `JEGKOW → nothing`, `CisPlatin → nothing` | measured on `main`, `OIN_EMIT_METAL_CONFIG=1` |

### The donor-SET finding, which is structural and worth its own row

**A distance heuristic cannot determine the donor set.** Metal-neighbour distance ratios
(`d / d_closest`), measured:

```
ZUMNEC (Mo):  O 1.00  1.00  1.02  1.02  1.03  1.03  |  C 1.39  1.39  1.41
JEGKOW (Rh):  C 1.00  P 1.22  N 1.23  I 1.50        |  O 1.63  C 1.72
```

ZUMNEC's donor/non-donor boundary demands a cutoff **below 1.39**; JEGKOW's iodide demands **at
least 1.50**. **No single value satisfies both.** Largest-relative-gap also fails — JEGKOW's
biggest gap (1.23 → 1.50) falls **before** the iodide, so it would drop a real donor. **Three
cutoffs were tried** (absolute 2.6 Å, absolute 3.0 Å, ratio 1.45) and each broke one fixture or
the other. **Resolution: take perception's metal-incident bonds from `get_tmc_mol`.**

### Refuted — three of the four formulations. See *Dead ends* for each.

### Two further bugs found during the lane, both worth recording

**(a) A test asserted a REGULAR tetrahedron is chiral.** It is not — Td contains improper
operations. And the old signed-volume descriptor **agreed with the bad test**, because a signed
volume is non-zero for *any* non-coplanar set of four **labelled** points. **"The labelling has
an orientation" is not "the shape is chiral"**, and only the permutation-invariant form can tell
them apart. Corrected: the regular tetrahedron is now an **achiral control**
(`test_a_REGULAR_tetrahedron_is_achiral`) with a **scalene** one as the chiral case
(`test_a_SCALENE_tetrahedron_is_chiral_and_inverts`).

**(b) `_admissible_permutations` initially yielded ZERO permutations for every input.**
`[itertools.permutations(x)] * n` repeats **one** iterator, so after the first is consumed the
rest are empty and `itertools.product` collapses to nothing. `is_achiral_chelate_aware`'s "no
symmetry found → chiral" verdict therefore came from a loop **that never ran** — it returned
"chiral" for *everything*, **including the achiral fixture, which is how it was caught**. The
first run was reported as "the first genuine detection of Δ/Λ"; it was **vacuous**. Fixed with
precomputed lists (`_INTERNALS`), with the trap named at the site, and pinned by
`TestAdmissiblePermutations::test_counts_are_exact` — because **an empty generator is
indistinguishable from a confident answer at the call site.**

---

## What was done

Nine commits on `main`, in order (all 2026-07-26), shipped in release `d799de1f`:

| commit | title |
|---|---|
| `26f504e3` | `lane5(v0.4.6): metal Delta/Lambda descriptor — sound, with two MEASURED blockers named` |
| `c72fbc1b` | `lane5(v0.4.6): SOLVE the homoleptic blocker — permutation-invariant chirality index` |
| `13f2b999` | `docs(lane5): the magnitude threshold is REFUTED — pucker looks like helicity` |
| `c172a57c` | `lane5: the descriptor's INPUT is wrong — donor positions cannot express Delta/Lambda` |
| `848af5c6` | `lane5: Delta/Lambda DETECTED — chelate-aware symmetry test, both fixtures correct` |
| `5b363057` | `test(lane5): put the chelate-aware Delta/Lambda descriptor under the suite` |
| `e8485ebf` | `lane5: WIRE the Delta/Lambda token to the emit path behind OIN_EMIT_METAL_CONFIG` |
| `820b92f5` | `lane5: teach the round-trip key to FOLD the \|mc:\| token, with the un-fold obligation recorded` |
| `27089512` | `lane5: generator now REPRODUCES the requested Delta/Lambda — Lane 5 complete end to end` |

### The notation decision: a trailing sidecar `|mc:+|` / `|mc:-|`, deliberately NOT `@OHn` / `@SPn`

**Why not the OpenSMILES non-tetrahedral tags:** they are defined against the **atom order in
the SMILES string**. That is *exactly* what makes RDKit's raw `_chiralPermutation`
non-reproducible across a re-parse, and it is the defect Lane 8 had to fix for tetrahedral tags.
A sidecar is self-describing instead.

Practical reasons the sidecar is cheap and safe, all following the landed ` |ax:±|` precedent:
it survives the parser (`generation/oin_parser.py` already strips sidecars), it moves no `{n}`
slot marker, it leaves `[El_GEO]` untouched so `METAL_REGEX` and `_METAL_STEREO_RE` are
unaffected, and it cannot collide with a ligand body.

### The descriptor — `src/oinsmiles/oin/metal_config.py`

| function | role |
|---|---|
| `chirality_index(donor_positions)` | the permutation-invariant Osipov–Pickup–Dunmur pseudoscalar, summed over **every** ordered 4-tuple and normalized by term count. Still used for the **sign** in the shipped path. O(n⁴) in donors; ≤ 8 donors ⇒ at most a few thousand terms |
| `metal_config_sign` / `metal_config_token` | formulation **2** (magnitude threshold, `_CHIRALITY_EPS = 1e-9`). **Retained but refuted as an achirality decision** |
| `_kabsch_proper_rmsd(a, b)` | RMSD after optimal **proper**-rotation superposition. Proper only — an improper "rotation" would map any set onto its own mirror and make every structure look achiral, which is the distinction being tested. SVD's last singular vector is flipped when `det < 0` |
| `is_achiral` | formulation **3** (unconstrained permutation search). **Retained, correct as a point-set test, refuted as a Δ/Λ test** |
| `metal_config_sign_symmetry` | formulation 3's token wrapper. Left in place, **unwired**, with the refutation written into its docstring |
| `_admissible_permutations(groups)` | formulation **4**'s core: relabellings that map whole chelates onto whole chelates of the same size |
| `is_achiral_chelate_aware(pts, groups, tol=0.35)` | formulation **4**. The form that works |
| `metal_config_token_chelate(pts, groups)` | the shipped token: `""` when chelate-aware-achiral, else `\|mc:+\|` / `\|mc:-\|` from the index's sign |
| `token_for_mol(mol)` | the encoder's one-line call site: derives **both** inputs from the mol — donors are the metal's perceived neighbours, the chelate partition is the connected components left after deleting the metal (which is what makes a bidentate's two donors *one* chelate). Returns `""` on anything unexpected rather than raising, because it runs inside the serialization path where an exception reroutes instead of surfacing |
| `parse_metal_config_token(oin)` | `\|mc:+\|` / `\|mc:-\|` / `None`. `None` for every OIN encoded without the lever — which is what makes the generator branch **inert by construction** |

Constants and why they are not tuned numbers:

- `_ACHIRAL_RMSD_TOL = 0.35` — chosen from a **24×** measured separation (1.3752 vs 0.0582 Å), so
  it is a wide band. Crystallographic pucker in an achiral sphere leaves a residual far below it,
  while a genuine Δ/Λ helix cannot be superimposed on its mirror by **any** chelate-preserving
  proper rotation.
- `_MAX_PERM_DONORS = 8` — above this, `is_achiral` returns `False` (**assume chiral**) rather
  than silently sampling, because a wrong "achiral" would suppress a real descriptor, which is
  the worse error. 8! = 40 320 cheap superpositions is affordable; 9! = 362 880 starts to matter
  in a hot encode path.
- `_CHIRALITY_EPS = 1e-9` — only has to clear floating-point residue now that the *decision* is
  made by symmetry rather than magnitude. Its first value, `1e-3`, was **larger than the real
  signal** (4.81e-4) and would have called ZUMNEC achiral.

### Emit — `src/oinsmiles/utils/xyz2mol.py`

Behind `OIN_EMIT_METAL_CONFIG`, **default OFF** (block at `:1390-1396`; registered in
`levers.py::_HELD_OFF`). Appended at the **single** return site, `:2038`:

```python
return inline_oin + _axial_suffix + _metal_config_suffix
```

**Computed BEFORE `_align_to_pai`**, beside the axial token and for the same documented reason:
principal-axis alignment **may reflect the coordinates**, and a reflection **inverts** a
chirality descriptor.

Note the correction recorded in `e8485ebf`: the author had earlier declined to wire it, claiming
it "touches the marker-placement path where a mistake silently corrupts coordination". That was
wrong — it is a *trailing sidecar*. **Overstating a risk is its own kind of inaccuracy.**

### Key fold — `src/oinsmiles/oin/compare.py`

`_METAL_CONFIG_TOKEN_RE = re.compile(r"\s*\|mc:[+\-]\|")`, applied in
`normalize_oin_for_comparison` alongside `_AXIAL_TOKEN_RE`. Folded **deliberately and
temporarily**: the token is a brand-new descriptor, so comparing it would turn every emitting
molecule into a round-trip failure the moment the lever is switched on — converting a silent
collapse into a loud one **without fixing anything**.

⚠ **The un-fold obligation is recorded at the regex itself**: the fold must be **removed in the
same commit** that promotes `OIN_EMIT_METAL_CONFIG`. **A key that folds an axis is not a valid
acceptance predicate for that axis** — and that is exactly why P1 was a blind spot in the first
place: *the key has always folded metal stereo, so no round-trip measurement could ever have
revealed the Δ/Λ collapse.* Adding a fold without recording its removal condition is how the
blind spot was created originally.

### Generator reproduction — `src/oinsmiles/generation/metallogen_adapter.py`

`target_mc = parse_metal_config_token(getattr(parsed, "original_oin", None))` computed beside
`target_axial` (~line 1489); the accept-first predicate then requires:

```python
target_mc is None or token_for_mol(cmol) == target_mc
```

**Why this check must exist separately from the key match** (~line 1554, in-code):

> `compare.py` FOLDS `|mc:|` (`_METAL_CONFIG_TOKEN_RE`), so a key match says nothing about
> helicity: accepting on the key alone would hand back the wrong enantiomer while reporting
> success.

**Inert by construction when the lever is off:** `parse_metal_config_token` returns `None` for
every OIN encoded without `OIN_EMIT_METAL_CONFIG` — the whole corpus by default — and `None`
short-circuits the condition, so selection stays byte-identical to pristine.

---

## Dead ends and refutations

### D1 — Signed volume of the first four donors in canonical slot order (formulation 1)

**What it was.** Order donors by a rotation/renumbering-independent key
(element, distance-to-metal, sorted multiset of distances to every other donor), take the signed
volume — normalized to a dimensionless triple product in [−1, 1] — and threshold its magnitude
at 0.15.

**It looked good.** ZUMNEC read **−0.862** against JEGKOW's **+0.026** — a **factor-30** margin,
so 0.15 was a wide band. It was invariant under 6 random proper rotations and inverted cleanly
under reflection.

**The measurement that killed it.** ZUMNEC is **homoleptic**: its six O donors are
symmetry-equivalent, so **every term in the ordering key is identical across them** — element,
distance to metal, and the multiset of distances to the other donors all tie. The ordering
therefore falls through to the order the atoms happened to arrive in, and **some of those
resolutions are related by an ODD permutation, which inverts a signed volume**. Measured:
**sign flipped 1 → −1 on 2 of 4 random renumberings.**

**Why no cleverer key fixes it.** For a homoleptic complex, *by the definition of
symmetry-equivalent*, **no scalar invariant separates the donors**.

**Sub-refutation recorded along the way:** an **absolute** volume threshold was wrong in *kind*,
not value — signed volume scales as (bond length)³, so no single number serves a 2.0 Å Rh–N
sphere and a 2.7 Å Rh–I one. That is why the dimensionless form exists at all.

**⚠ Related trap, same family:** this is the *same* trap that made the Y2 axial token
reflection-invariant — canonicalizing symmetry-equivalent elements by a value that the reflection
also changes.

### D2 — Magnitude threshold on a permutation-invariant chirality index (formulation 2)

**What it was.** Stop needing an ordering: sum the Osipov–Pickup–Dunmur pseudoscalar over
**every** ordered 4-tuple. Relabelling permutes the terms without changing the total.

**The parts that are real and remain real** (measured on ZUMNEC, the fixture the ordered version
failed): 6 random donor permutations give an **identical** value; 3 proper rotations identical;
reflection gives **exact negation**; a perfect square and an ideal octahedron both return
**exactly `+0.000e+00`** — not rounding noise.

**The measurement that killed it.** The accompanying claim — "achirality falls out of the index,
no planarity test needed" — was generalized from **two synthetic controls**. Re-measured with
donors taken from **perception**:

| fixture | index |
|---|---|
| ZUMNEC — genuinely chiral Δ/Λ | **−4.807e-04** |
| JEGKOW — achiral square planar, only puckered | **−3.287e-04** |
| ideal square (synthetic) | +0.000e+00 |

**The two real structures are the same order of magnitude, 1.5× apart.** Crystallographic pucker
in an *achiral* complex produces a chirality index comparable to genuine helicity, so **no
threshold separates them** and any `_CHIRALITY_EPS` is arbitrary whatever value it takes. The
exact-zero cancellation is a property of **idealized coordinates only**.

**Second error the same commit exposed:** `_CHIRALITY_EPS` was first set to **1e-3** — *larger
than the real signal* (4.81e-4) — so it would have called ZUMNEC achiral and silently emitted
nothing for the one fixture the lane exists for.

**The transferable lesson, in the author's own words:** *a measurement that only exercised the
easy case (synthetic ideal geometry) confirmed a wrong belief. The fixtures were available the
whole time; I read the controls first and generalized.*

### D3 — Unconstrained mirror-superposition symmetry test (formulation 3)

**What it was.** Chirality is a **symmetry** property, so replace the magnitude threshold with a
symmetry test: mirror the donor set, then ask whether **any** relabelling of the mirrored points
can be superimposed on the original by a **proper** rotation. If one can, the mirror is the same
object.

**The measurement that killed it.**

| fixture | symmetry test says | truth |
|---|---|---|
| ZUMNEC — chiral Δ/Λ tris-bidentate | **achiral** | chiral |
| JEGKOW — square planar | achiral | achiral |
| ideal square | achiral | achiral |

**And it is not the test that is wrong — it is the input.** As a *bare point set* ZUMNEC **IS**
achiral: six oxygens at octahedral vertices admit improper operations; there is no handedness in
donor positions alone. A permutation search over **unlabelled** points is free to **re-pair**
donors into different ligands, **which no physical operation can do**, so it always finds a
"symmetry" that is not chemically available.

**It also retro-explains D2's apparent success:** `chirality_index`'s non-zero reading for ZUMNEC
(−4.807e-04) was **residual crystallographic distortion**, the same magnitude as achiral
JEGKOW's pucker (−3.287e-04). **It was never detecting helicity at all.** Three of its four
proven properties (rotation invariance, reflection inversion, permutation invariance) were real
and remain real — they were properties of a **quantity that does not mean what the lane needs**.

This is the only refutation in the release that invalidated the **input** rather than the method.

### D4 — Distance-cutoff donor-set determination

**Killed by** the ratio table above: ZUMNEC demands a cutoff **below 1.39**, JEGKOW's iodide
demands **at least 1.50**; largest-relative-gap fails because JEGKOW's biggest gap falls
*before* the iodide. Three cutoffs tried (2.6 Å, 3.0 Å, ratio 1.45); each broke one fixture.

**Residue in the tree:** `tests/unit/test_metal_config.py` still carries a distance-based donor
finder (`_donors_in_canonical_order`, `_DONOR_RATIO = 1.20`) tuned for ZUMNEC only, and
`TestJegkowSquarePlanarEmitsNothing::test_emits_nothing` is marked `@unittest.expectedFailure`
**as a harness limitation, not a descriptor defect** — the docstring records that with the right
donor set (C, P, N, I) the normalized triple product is **+0.026**, comfortably inside the
planarity band. The shipped path (`token_for_mol` → `_donors_and_chelate_groups`) uses perception
and is correct.

### D5 — "A regular tetrahedron is chiral"

A test asserted it; it is false (Td contains improper operations). The old signed-volume
descriptor **agreed with the bad test** because a signed volume is non-zero for any non-coplanar
set of four **labelled** points. Confusing "the labelling has an orientation" with "the shape is
chiral" is precisely the flaw the permutation-invariant form removes.

### D6 — A permutation generator that yielded nothing (the vacuous detection)

`[itertools.permutations(x)] * n` repeats **one** iterator ⇒ `itertools.product` collapses ⇒
`_admissible_permutations` yielded **zero** permutations for every input ⇒
`is_achiral_chelate_aware`'s "no symmetry found → chiral" came from a loop that never ran. Every
fixture came back chiral **including the achiral one, which is how it was caught**. The first run
was written up as "the first genuine detection of Δ/Λ" and **retracted**.

Same family as the empty-corpus and buffered-stdout failures already recorded in this project's
docs: **nothing measured, confident answer.**

---

## Where it landed

**Complete pipeline, end to end:** descriptor → emit (default-OFF lever) → key fold (with the
un-fold obligation recorded) → generator reproduction.

| stage | location | state |
|---|---|---|
| descriptor | `src/oinsmiles/oin/metal_config.py` | chelate-aware symmetry test (formulation 4) |
| donors + chelate partition | `token_for_mol` via `get_tmc_mol` perception | from perception, not coordinates |
| emit | `src/oinsmiles/utils/xyz2mol.py:1390-1396`, appended at `:2038` | behind **`OIN_EMIT_METAL_CONFIG`**, **default OFF** |
| key fold | `src/oinsmiles/oin/compare.py::_METAL_CONFIG_TOKEN_RE` | folds; **must un-fold on promotion** |
| generator | `src/oinsmiles/generation/metallogen_adapter.py` (`target_mc`, accept-first) | inert when `target_mc is None` |

**Verified emit behaviour** (measured for this retrospective on `main` with
`OIN_EMIT_METAL_CONFIG=1`):

```
ZUMNEC              ... 1O{5} |mc:-|
ZUMNEC (z-mirror)   ... 1O{4} |mc:+|      <- INVERTS, as a chirality descriptor must
JEGKOW              (no mc token)         <- square planar, achiral: correctly silent
CisPlatin           (no mc token)         <- achiral control: no false positive
```

**Guard tests — `tests/unit/test_metal_config.py`, 16 tests, `OK (expected failures=1)`:**

| class | tests | what it pins |
|---|---:|---|
| `TestZumnecHelicity` | 4 | non-planar & emits; invariant under 6 proper rotations; **inverts under reflection** (the property Y2 failed); invariant under 4 atom renumberings (*was* an xfail — solved by dropping the need for an ordering) |
| `TestJegkowSquarePlanarEmitsNothing` | 1 (**xfail**) | harness donor-finder limitation, named rather than papered over |
| `TestDegenerateInputs` | 4 | `<4` donors silent; exactly coplanar silent; **REGULAR tetrahedron achiral**; **SCALENE tetrahedron chiral and inverts** |
| `TestAdmissiblePermutations` | 3 | `test_counts_are_exact` (**48 / 4 / 24**) — the generator is **not empty**; every yield is a genuine permutation; chelate membership preserved |
| `TestChelateAwareDeltaLambda` | 4 | ZUMNEC chiral and the token **inverts**; JEGKOW emits nothing; **`test_unconstrained_search_disagrees_on_zumnec`** (pins *why* the partition is required, not merely that it helps); invariant under 4 proper rotations |

**Fixture guards — `tests/unit/test_metal_stereo_fixtures.py`, 11 tests,
`OK (expected failures=1)`:** geometry cleanliness gate (0 vdW clashes, single metal centre);
`TestZumnecTrisBidentate` (three chelates with symmetry-equivalent donors ⇒ no fac/mer to lean
on; mirror distinct > 1.0 Å; `@OH` permutation flips; no other stereo axis implicated;
**`test_currently_key_blind`** asserts the P1 collapse still reproduces on the default path);
`TestZumnecAspirational::test_metal_chirality_should_diverge_at_key` (**xfail** — flips to an
unexpected success when the lever is promoted and the fold removed);
`TestJegkowSquarePlanar::test_mirror_is_the_same_isomer` (correct invariance, not a weak fixture).

**Also still xfail on the default path:**
`tests/unit/test_injectivity_probes.py::test_metal_chirality_should_diverge_at_key`.

### ⚠ What is NOT usable in the shipped default, and why

**On the shipped default configuration the encoder emits no metal-configuration descriptor at
all, so P1 remains exactly as blind as Y1 measured it.** Three separate reasons, all deliberate:

1. **`OIN_EMIT_METAL_CONFIG` is default OFF** — the standard information-ADDING trade. The
   generator must be able to reproduce what the encoder emits; promoting converts a *silent*
   collapse of Δ/Λ enantiomers into a *loud* round-trip failure. Promotion is an explicit product
   call, not something to take inside an encoder change.
2. **The round-trip key still FOLDS the token.** So even with the lever on, the batch harness
   cannot confirm the descriptor — only the raw string can. Conversely, promoting the lever
   without removing the fold would leave the round trip **structurally unable to verify the one
   thing the token encodes**.
3. **No corpus population measurement exists.** 0/150 was the *pre-descriptor* count. Nobody has
   measured how many corpus molecules would emit a `|mc:±|` token.

### ⚠ Stale docstrings inside `metal_config.py` — read the tests, not the module header

The module was written incrementally and **three of its docstring claims were later refuted by
its own author but not deleted**:

| location | stale claim | status |
|---|---|---|
| module docstring, "THE DESCRIPTOR" bullet 4 | *"exactly 0 for an achiral arrangement — a perfect square and an ideal octahedron both return `+0.000e+00`. Achirality falls OUT of the index instead of needing a planarity test beside it, which is why square-planar `JEGKOW` emits nothing without a special case."* | **REFUTED** by D2 — true only for idealized coordinates; real JEGKOW reads −3.287e-04 |
| `_CHIRALITY_EPS` comment | *"an achiral point set cancels to **exactly 0.0** … The separation is therefore between exact zero and the signal"* | **REFUTED** by D2 |
| module docstring, "STATUS — DESCRIPTOR AND VALIDATION ONLY, NOT WIRED TO EMIT" | *"deliberately not yet called from `xyz2mol.py`'s emit path, and no lever turns it on"* | **STALE** — `e8485ebf` wired it; `xyz2mol.py:1391` calls `token_for_mol` behind `OIN_EMIT_METAL_CONFIG` |

The correct, current account lives in `is_achiral`'s and `metal_config_sign_symmetry`'s
docstrings, in `_ACHIRAL_RMSD_TOL`'s comment, and in `docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md` (the three
"Lane 5, Nth measurement" sections). **Do not quote the module header.**

---

## Open questions / for the next agent

### The promotion gate for `OIN_EMIT_METAL_CONFIG` — three conditions

1. **Measure the corpus population.** Nobody has run `token_for_mol` over the dataset. The
   *only* number in existence is the pre-descriptor **0/150** metal-`@` count, which says nothing
   about how many molecules have ≥4 donors, a chelate partition, and a chelate-aware-chiral
   sphere. This is the cheapest missing measurement in the lane and it is **load-independent**
   (string/geometry only, no wall-clock), so it can be taken while a sweep runs. Model it on
   `tools/injectivity/axial_population.py`, which produced the equivalent 2.49% figure for P2.
2. **Remove the key fold in the same commit.** `compare.py::_METAL_CONFIG_TOKEN_RE` must stop
   folding when the lever is promoted. Not before: changing it while the lever is off is inert,
   and changing it with the lever off but a token present would silently tighten the key for
   anyone experimenting. The two changes belong in **one** commit.
3. **Demonstrate generator reproduction as a rate, not a mechanism.** The accept-first hook
   exists and is correct, but there is **no A/B** for helicity comparable to the axial lane's
   22/22-vs-8/22. Build one on the `axial_cohort_ab.py` pattern: encode both enantiomers with the
   lever ON, generate from each OIN, read the generated structure's own `token_for_mol`, and
   report a match rate for a selection arm against a no-selection arm. **Two fixtures is the
   sample size that produced four wrong answers elsewhere in this release.**

### Specific open technical questions

- **Does `_ACHIRAL_RMSD_TOL = 0.35` survive a corpus?** It is justified by a 24× separation on
  **two** structures. The tolerance sits between crystallographic pucker in an achiral sphere and
  a genuine helix, and nothing has probed the middle of that band. A histogram of best
  chelate-preserving mirror-superposition RMSD over the corpus would show whether the gap is a
  real bimodality or an n=2 artefact.
- **What happens above 8 donors?** `is_achiral` returns `False` (assume chiral) above
  `_MAX_PERM_DONORS = 8`, deliberately choosing the safer error. But
  `is_achiral_chelate_aware` has **no such guard** — it enumerates
  `_admissible_permutations(groups)` unconditionally, and `_INTERNALS` is only precomputed for
  group sizes **1–6**. A chelate of denticity 7+ would raise a `KeyError` inside the generator,
  which `token_for_mol`'s blanket `except Exception: return ""` would swallow into silence. Worth
  either a guard or an explicit test.
- **Square-planar diastereomers are still uncovered.** JEGKOW correctly emits nothing, because
  reflection is the wrong distinctness operator for `@SP` — the coordination plane *is* a mirror
  plane. The right operator is a **donor swap**, which permutes which donors sit *trans*. Lane 7
  built `swap_donor` in `tools/injectivity/twin_operators.py` and an octahedral donor-swap probe
  (`b315b929`), but **no `@SP` diastereomer descriptor exists**. That is a separate, unopened
  blind spot, not part of P1.
- **Clean up the refuted docstrings** listed above before someone builds on them. Three separate
  claims in `metal_config.py` assert properties measurement destroyed.
- **Two dead formulations are still exported.** `metal_config_sign`, `metal_config_token` and
  `metal_config_sign_symmetry` remain in `__all__`. They are deliberately retained (the
  pseudoscalar's three invariances are real and reused; the point-set achirality test is correct
  as such), but a future caller could pick the wrong one. Their docstrings carry the warnings;
  consider a stronger signal.
