# Lane: valence order (`swimlane/v045-lane2`, `OIN_STABLE_METAL_AC`)

**What the lane was for:** the adjacency-matrix perception in `xyz2AC_obabel` capped
over-valent atoms by iterating **input atom index order**, so the molecule the encoder
perceived depended on the order the atoms happened to be listed in the XYZ file. This lane
found the line, replaced the iteration order with a deterministic one derived from the
structure, and gated it — then it was promoted to default-ON.

> ## ⚠ Source discrepancy, flagged explicitly
>
> The retrospective brief lists `docs/agentic-notes/v0.4.5/VALENCE_ORDER_v0.4.5.md`, `tests/unit/test_valence_order.py`,
> `tools/valorder_probe.py`, `tools/valorder_feasibility.py`, `tools/valorder_encode_ab.py`
> and `git log main..swimlane/v045-valorder` as the sources for this lane. **They document a
> different thing.** Those artefacts are about the order in which `AC2BO` enumerates
> **candidate valence assignments** inside the bond-order search — a downstream, and purely
> internal, enumeration order. Their content is covered in
> `docs/agentic-notes/v0.4.5-retrospective/LANE-valence-search.md`, where it belongs.
>
> The lane described *here* — the atom-index order of the **valence-capping loop**, `DUDREA_comp_0`,
> the `OIN_STABLE_METAL_AC` lever — lives on `swimlane/v045-lane2` and is documented in
> `docs/agentic-notes/v0.4.5/RENUMBERING_INSTABILITY_v0.4.5.md` (§"Follow-up, measured 2026-07-25"),
> `docs/agentic-notes/v0.4.5/PROMOTION_GATE_v0.4.5.md`, `docs/agentic-notes/v0.4.5/CANONICAL_OIN_v0.4.5.md` §"Known gaps" item 5,
> commit `8bf9df61`, `src/oinsmiles/utils/xyz2mol_local.py:1858-1913`,
> `tools/geometry_tag_shift.py`, and `tests/unit/test_xyz2mol_errors.py`. Those are the
> sources this file is built from. Two lanes, two different meanings of the word "order" —
> read both files together and do not assume one is a stale copy of the other.
>
> Also note: `git log --oneline main..swimlane/v045-lane2` is **empty**. The branch is an
> ancestor of `main`; its work is released.

---

## ELI5

An XYZ file is just a list of atoms and their positions — and the list has an order, which is
arbitrary. The encoder decides which atoms are bonded by comparing distances (order doesn't
matter for that), but then it has to fix up atoms that ended up with too many bonds, and it did
that fix-up by walking the atom list from top to bottom, deleting the longest bond of whichever
over-bonded atom it met first. That means the outcome depended on the list order: **take the
same molecule, shuffle the lines in the file, and the encoder perceives a genuinely different
molecule.** On `DUDREA_comp_0` — a yttrium borohydride where one hydrogen bridges between the
boron and the metal — trimming the metal first leaves the Y–H bond intact and the metal reads as
5-coordinate square-pyramidal (`[Y_SPY]`); trimming that hydrogen first deletes the Y–H bond and
the same metal reads as 4-coordinate tetrahedral (`[Y_TET]`). Same coordinates, same file,
different answer, purely because of line order. The fix makes the trim order a property of the
structure (heaviest atom first, then invariant tie-breaks) instead of a property of the file.

---

## The work, visually

```
tests/.../DUDREA_comp_0.xyz            SAME COORDINATES, atoms listed in two orders
  ORDER A                                ORDER B
  0: Y   (metal)                         0: H   (the bridging hydride)
  1: H   (bridging hydride)              1: Y   (metal)
  2: B                                   2: B
  ...                                    ...

                         |                              |
                         v                              v
=====================================================================================
STEP 1 — the DISTANCE pass                       xyz2mol_local.py:1848-1856
         (inside xyz2AC_obabel, :1809)
=====================================================================================
    for i in range(N):
      for j in range(i+1, N):
        if dMat[i,j] <= Rcov_i + Rcov_j + tolerance:   AC[i,j] = AC[j,i] = 1

    A SYMMETRIC comparison of the distance matrix against covalent-radius sums.
    ORDER-FREE: permuting the atoms permutes the matrix and nothing else.

    Result for DUDREA, in BOTH orders (identical up to relabelling):
                                B
                               /|\        the bridging H is bonded to BOTH B and Y
                          H---H-H         -> H has degree 2, and H's max valence is 1
                          |    |          -> SOMETHING must be deleted
                          Y----H
       Y has six hydride contacts: 2.298 / 2.300 / 2.328 / 2.379 / 2.408 / 2.421 A

=====================================================================================
STEP 2 — the VALENCE-CAPPING loop                xyz2mol_local.py:1889-1913
         *** THE ONLY ORDER-DEPENDENT STEP IN AC PERCEPTION ***
=====================================================================================

  BEFORE (default pre-v0.4.5):        cap_order = range(num_atoms)   <-- INPUT ORDER
  AFTER  (OIN_STABLE_METAL_AC):       cap_order = sorted(range(num_atoms), key=_cap_key)

    for i in cap_order:
        if i in exempt: continue            <-- [boron lane composes HERE, see below]
        while sum(AC[i,:]) > max(atomic_valence[Z_i]):
            AC = remove_weakest_bond(mol, i, AC, dMat, pt)   # deletes i's LONGEST bond

  Why order decides the answer: capping atom i DELETES a bond, which LOWERS some other
  atom j's count -- so whether j still needs capping at all depends on whether i came first.

   ORDER A: Y visited first            |   ORDER B: the bridging H visited first
   ------------------------------------|--------------------------------------------
   Y is over-valent -> trim Y's        |   H is over-valent (degree 2 > 1)
   longest contacts first.             |   -> trim H's longest bond = the Y-H bond
   The Y-H bond SURVIVES.              |   The Y-H bond is DESTROYED.
   H's count drops to 1 -> H no        |   Y's count drops -> Y needs less capping.
   longer needs capping at all.        |
            |                          |            |
            v                          |            v
   metal coordination degree 5         |   metal coordination degree 4
   template fit -> SPY                 |   template fit -> TET
   emitted:  [Y_SPY]. ... .[BH3]{2}    |   emitted:  [Y_TET]. ... .[BH3]{2}
             .[H]{3}.[BH3]{4}.[H]      |             .[BH3]{3}.[H].[H]
   ------------------------------------|--------------------------------------------
       ^ A DIFFERENT PERCEIVED MOLECULE FROM THE SAME COORDINATES.
         The [M_XXX] tag is not cosmetic: it selects the vertex table -> the rotation
         group -> the canonical slot labelling -> the comparison key's whole signature.

   NOTE: the losing contact was NOT the shortest. It was the 2.328 A hydride, third of
   six. So this was never "shortest wins" -- it was iteration order.

=====================================================================================
THE FIX — where it intervenes                    xyz2mol_local.py:1894-1903
=====================================================================================
   if _lever_enabled("OIN_STABLE_METAL_AC"):
       def _cap_key(i):
           nbrs  = np.nonzero(AC[i, :])[0]
           dists = tuple(sorted(round(float(dMat[i, j]), 4) for j in nbrs))
           return (-atomic_number(i),   # (1) HEAVIEST ELEMENT FIRST
                   -len(nbrs),          # (2) then highest degree
                   dists)               # (3) then sorted neighbour distances
       cap_order = sorted(range(num_atoms), key=_cap_key)
   else:
       cap_order = range(num_atoms)                       # historical path, byte-identical

   Every component of the key is a ROTATION- and PERMUTATION-INVARIANT scalar, so
   genuinely equivalent atoms tie and everything else separates -- with no reference
   to an atom index anywhere.

   Effect on DUDREA: Y (Z=39) sorts before every H (Z=1) in EVERY input order, so
   ORDER A's outcome is always taken.  Measured: AC differed in 3/8 random renumberings
   with the lever off, 0/8 with it on.

=====================================================================================
COMPOSITION — two independent defects live in this ONE loop, and v0.4.5 fixed both
=====================================================================================
   exempt   = boron_cage_vertices(atomic_nums, AC)     if OIN_BORON_CAGE   [WHICH atoms]
   cap_order = sorted(..., key=_cap_key)               if OIN_STABLE_METAL_AC  [what ORDER]
                     |
                     v
   for i in cap_order:  if i in exempt: continue   ...
   ^^^ They compose in the same statement and must be read together.

LEGEND
  AC ................... adjacency matrix; AC[i,j]=1 means i and j are bonded
  dMat ................. 3D distance matrix (Chem.Get3DDistanceMatrix)
  Rcov ................. covalent radius; tolerance defaults to 0.45 A
  cap / capping ........ deleting an atom's longest bonds until its degree is within
                         max(atomic_valence[Z])
  [M_XXX] .............. the metal geometry token in the OIN string, e.g. [Y_SPY]
  {n} .................. an OIN coordination slot marker
  exempt ............... vertices the boron-cage lane withholds from capping
  _cap_key ............. the invariant sort key this lane introduced
```

---

## Initial assumptions and hypothesis

The lane did not start from this defect. It started from a **measurement nobody had taken**:
`tools/canonicality_probe.py` re-presents one input structure three ways while holding the
molecular graph fixed — `rotate` (a random **proper** rotation, det = +1), `renumber` (the order
atoms appear in the XYZ file), and `both` — so the correct answer is byte-**identical**, a known
ground truth rather than an inferred baseline. Rotations are forced proper because an improper
operation mirrors the structure, which legitimately changes a chiral molecule's encoding.

Measured at `main` @ `20044883`, all levers OFF, 225 molecules (seed 42) from the
25,197-basename corpus, 2 trials per transform, 223 encoding successfully:

| | count | share |
|---|---:|---:|
| byte-stable across all transforms | 125 / 223 | **56.1%** |
| drifted | 98 / 223 | **43.9%** |
| **of which the comparison KEY also changed** | 47 / 223 | **21.1%** |

Drift by transform: `renumber` 132, `both` 123, **`rotate` 0**. Severity classes under pure
renumbering: stereo-only flip **29 (13.0%)**, other drift 79 (35.4%), **geometry
classification changed 1 (0.4%)**, aromaticity perception collapsed (>2 atoms) 8 (3.6%).

The going-in hypotheses, in the order the probe's own scope call listed them:

1. `AC2BO` / `get_UA_pairs` order-dependence changing perceived bond orders, hence CIP
   priorities — believed to share a root with Lane 1's `OIN_CANONICAL_PERCEPTION`;
2. `core/chirality.py`'s `CIPAssigner` / `ChiralityRecoveryUtility` ordering assumptions;
3. `_align_to_pai`'s index-dependent pivot and the `(i+1)**3` Z-moment sign
   (`utils/xyz2mol.py:941`, `:971`) — nominated as the most likely route for `DUDREA_comp_0`,
   on the reasoning that the Z-moment weighting flips Y/Z under renumbering and the geometric
   template fit then selects a different polyhedron.

Lane 2 also arrived carrying a design of its own: a **canonical slot post-pass**, on the premise
that geometry-tag drift was a slot-labelling problem.

Three of those four beliefs were wrong about `DUDREA`.

---

## What was actually found

### CONFIRMED — the encoder is fully orientation-invariant

Across 225 molecules × 2 trials, **not one** structure changed its OIN string under a random
proper rotation. `rotate` drift is **0** in every arm of every subsequent measurement, including
the final promotion gate. `_align_to_pai` does its job for orientation. **Every defect in this
programme is an atom-numbering dependence, not a geometry one** — which is what made the
`(i+1)**3` Z-moment hypothesis worth dropping.

### CONFIRMED — the defect is one loop, and it is not the one anybody named

`xyz2AC_obabel`'s distance pass is order-**free**: a symmetric comparison of the distance matrix
against the covalent-radius sums (`src/oinsmiles/utils/xyz2mol_local.py:1848-1856`; the lane
docs cite the pre-merge location `:1194-1202`). **The only order-dependent step in AC perception
is the valence-capping loop that follows** (`:1889-1913`), which iterated `for i in
range(num_atoms)` — in input atom order.

The mechanism, stated exactly: capping atom *i* removes a bond, which lowers some atom *j*'s
count, so whether *j* still needs capping **depends on whether *i* was visited first**.

### CONFIRMED — the cost, on the clearest available example

`DUDREA_comp_0` is a Y borohydride whose bridging hydride is bonded to **both** B and Y,
exceeding H's valence of 1, so something must be deleted:

* cap **Y** first → the Y–H bond survives → metal coordination degree **5** → geometry **SPY**;
* cap that **H** first → it drops the Y–H bond instead → degree **4** → geometry **TET**.

Emitted strings, from the probe:

```
BASE : [Y_SPY]. ... .[BH3]{2}.[H]{3}.[BH3]{4}.[H]
RENUM: [Y_TET]. ... .[BH3]{2}.[BH3]{3}.[H].[H]
```

Not only the tag changed — **a different set of atoms is recorded as coordinated.** That is an
isomer-level change produced by nothing but atom numbering, and it moves the comparison key.

**The losing contact was not even the shortest.** The metal has six hydrides at 2.298 / 2.300 /
2.328 / 2.379 / 2.408 / 2.421 Å, and the one that flipped was **2.328 Å** — third of six. So
this was never "shortest wins"; it was iteration order.

**Measured:** the perceived adjacency matrix differed in **3 of 8** random renumberings with the
fix off, and **0 of 8** with it on.

> ⚠ **A numeric inconsistency in the sources, recorded rather than smoothed over.** The
> SPY/TET narrative is stated in terms of metal coordination degree **5 vs 4** (docs and commit
> body both), while the same commit's measurement paragraph says the metal degree was **flipping
> between 6 and 7** off the lever and settles **consistently at 7** with it on. The most likely
> reading is that 5/4 describes the coordination count the geometry template fit consumes
> (distinct `{n}` donor slots) while 6/7 is the raw AC row sum on the metal, but **this was not
> verified**. If you need the exact relationship, re-derive it; do not quote both numbers as if
> they measure the same quantity.

### REFUTED — that Lane 1's perception lever would close the stereo class

This hypothesis had been reasoned out and **communicated in writing to the Lane 8 agent**: the
13% stereo-flip class was assumed downstream of order-dependent bond-order perception and
therefore assumed to be fixed by `OIN_CANONICAL_PERCEPTION`. Measured on the three worked
examples, 3 trials × 3 transforms:

| levers | byte-stable | key-level defects remaining |
|---|---|---|
| all OFF | 0/3 | **3** |
| `OIN_STABLE_METAL_AC` | 0/3 | **2** |
| `+ OIN_CANONICAL_PERCEPTION` | 0/3 | 2 |
| `+ OIN_CANONICAL_BODY + OIN_CANONICAL_SLOTS` (all four) | 0/3 | **2** |

**Only the metal-AC fix closes anything.** `FEQFIS_comp_0`'s stereo flip and `CEBVIR_comp_0`'s
aromaticity collapse survive all four levers. So Lane 8's independent investigation was
genuinely required rather than redundant. (Scope caveat carried in the source: n = 3,
deliberately the hardest hand-picked cases; it refutes the *specific* claim, not Lane 1's
general result of 6 fixed / 0 regressed over 250 molecules.)

The stereo class was subsequently closed by Lane 8's own `OIN_STABLE_STEREO` — over 10 of the 29
known stereo-flip molecules, byte-stable went **0/10 → 8/10** and key-level defects **10 → 1** —
with a mirror guard confirming 10/10 mirrors still produce a different string (`mirror ==
swap(base)`, `@`↔`@@`), so nothing collapsed. That is a sibling lane's result, recorded here only
because this lane's refutation is what made it necessary.

### REFUTED — that this was a slot-labelling problem

`DUDREA_comp_0` "was never a slot-labelling problem, which is why Lane 2's slot post-pass could
not touch it." The lane's own planned design was refuted by its own diagnosis.

### CONFIRMED — the corpus is clean under the change (the veto that was set against the fix)

The direction of the change is asymmetric and dangerous: **capping the metal first can only ADD
bonds the old order discarded.** More metal bonds means a higher coordination number, which
means the geometric template fit can land on a *different polyhedron* — corpus-wide, silently.
And `[M_XXX]` is not cosmetic: it selects the vertex table, hence the rotation group, hence the
canonical slot labelling, hence the comparison key's entire vertex signature. Fixing one
molecule cleanly is not evidence of safety at scale.

`tools/geometry_tag_shift.py` was written specifically to enforce that veto. Per molecule it
runs one encode with the lever off and one with it on, then reports the `[M_XXX]` transition
matrix (off-tag → on-tag, so a shift is visible **by direction**, not only as a count),
coordination number taken as the number of distinct `{n}` slots (so a tag change can be
attributed to a genuine donor gain rather than a template tie flipping), and whether the change
is confined to the tag or the whole string moved.

**Result — PASS, 298 molecules: 0 string changes, 0 `[M_XXX]` changes, 0 coordination-number
changes, no transitions.** `0/298` refutes the concern the veto existed for.

### CONFIRMED — the promotion evidence

| measurement | result |
|---|---|
| all six canonicality levers, byte-stability under rotation/renumbering (300-molecule seed-42 sample, 298/299 encoded) | **58.1% → 69.6%** (+11.5 pts, +35 molecules) |
| comparison-key instability | **60 → 16 molecules** (1 in 5 → roughly 1 in 19; a 73% reduction) |
| re-baseline over 936 molecules | **145 of 436** previously-failing molecules **FIXED (33.3%)**; of 500 previously-passing guards, all **11** apparent regressions are `TimeoutException exceeded 300s` against a capstone baseline that ran at 1800 s ⇒ **zero correctness regressions** |
| `rotate` drift | **0 in both arms** — orientation-invariance preserved |
| goldens on the opt-out path | byte-identical (`test_regression_stability.py`) |
| suite | **837 tests OK** (3 skipped, 4 expected failures) at release |

> ⚠ **Attribution caution.** The 145-fixed re-baseline and the 58.1% → 69.6% figures were
> measured with **all six** promoted levers ON, not with `OIN_STABLE_METAL_AC` alone. The
> docstring of `tests/unit/test_xyz2mol_errors.py` cites "capstone A/B: 145 molecules fixed,
> zero correctness regressions" in support of this lever specifically; treat that as *the set
> containing this lever showed no correctness regressions*, which is what was measured. The
> lever-specific evidence is `geometry_tag_shift` **0/298** and DUDREA's **3/8 → 0/8**.

---

## What was done

### The fix

`src/oinsmiles/utils/xyz2mol_local.py`, in `xyz2AC_obabel` (`:1809`), 41 insertions / 4
deletions (commit `8bf9df61`). The index order is replaced by a canonical one:

```python
if _lever_enabled("OIN_STABLE_METAL_AC"):

    def _cap_key(i):
        nbrs = np.nonzero(AC[i, :])[0]
        dists = tuple(sorted(round(float(dMat[i, j]), 4) for j in nbrs))
        return (-mol.GetAtomWithIdx(int(i)).GetAtomicNum(), -len(nbrs), dists)

    cap_order = sorted(range(num_atoms), key=_cap_key)
else:
    cap_order = range(num_atoms)
```

Key design points, each with its reason:

* **`-atomic_number` — heaviest element first.** Chemically the right way round: a κ²/κ³ BH₄
  really *is* bound through its hydrides, so the metal should claim them before hydrogen's
  valence rule discards them. It also happens to reproduce **today's** answer for `DUDREA`,
  where the metal was atom 0 in the shipped input — so the fix keeps the current output for
  that molecule rather than changing it and stabilising the change.
* **`-len(nbrs)` then `sorted(neighbour distances)` — invariant tie-breaks only.** Both are
  rotation- **and** permutation-invariant scalars, so genuinely equivalent atoms tie and
  everything else separates **without any reference to an atom index**. Distances are rounded to
  4 decimal places so floating-point noise cannot reorder equivalent atoms.
* **No index fallback.** A final `i` tie-break would have reintroduced exactly the dependence
  being removed; ties are left as ties.

### Gating and promotion

* Lever name: **`OIN_STABLE_METAL_AC`**. Shipped **default OFF** in `8bf9df61` — deliberately,
  because the change is to **perception**, not serialization, and *"do not promote on the
  strength of one fixture."*
* Promoted to **default-ON** in the v0.4.5 Wave D promotion, registered in
  `src/oinsmiles/oin/levers.py::_DEFAULT_ON` alongside `OIN_CANONICAL_BODY`,
  `OIN_CANONICAL_PERCEPTION`, `OIN_CANONICAL_SLOTS`, `OIN_CANONICAL_ETA_WINDING` and
  `OIN_STABLE_STEREO` (`OIN_BORON_CAGE` joined the set in v0.4.6).
* Read through `lever_enabled(name, override=None)`, where `"0"`, `"false"`, `"no"`, `"off"` and
  `""` disable. This matters: the older bare-truthiness spelling `os.environ.get("X")` made
  `X=0` **enable** X, and `OIN_BORON_CAGE` alone had five sites on that spelling.
* The promotion pattern is the `OIN_EARLY_EXIT` template from v0.4.4
  (`metallogen_adapter.py:1636-1653`): membership test on `ff_params` so an explicit `False` can
  opt out, and an in-code comment naming the evidence document.

### Why it was safe to promote *with* the other five

The rule the release adopted, and the one worth keeping: **each of the six repairs a renumbered
presentation without rewriting the canonical answer.** That is why the corpus shows no churn.
Levers that **add information** to the string (`OIN_EMIT_AXIAL`, `OIN_EMIT_LOCKED_DONOR`,
`OIN_EMIT_METAL_CONFIG`) are a different trade — the generator must then be able to reproduce
what they emit, so promoting one converts a silent false positive into a loud false negative.
Those stayed opt-in, with reasons recorded in `levers.py::_HELD_OFF`.

### Composition with the boron-cage lane — read the two together

Both defects live in the **same** capping loop and v0.4.5 fixed both, so the statement is now
two-part:

```python
exempt = set()
if _lever_enabled("OIN_BORON_CAGE"):
    atomic_nums = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in range(num_atoms)]
    exempt = boron_cage_vertices(atomic_nums, AC)     # WHICH atoms get capped

if _lever_enabled("OIN_STABLE_METAL_AC"):
    ...
    cap_order = sorted(range(num_atoms), key=_cap_key)  # what ORDER they are capped in
else:
    cap_order = range(num_atoms)

for i in cap_order:
    if i in exempt:
        continue
    ...
```

`OIN_BORON_CAGE` exists because the loop deletes bonds while connectivity exceeds
`max(atomic_valence[Z])`, and for boron that cap is 4 while a closo/nido deltahedral vertex has
5–6 neighbours — so on a carborane the loop amputates 7–19 B–B cage edges, shattering an intact,
correctly-perceived cage. The exemption is scoped to element B in a B–B–B triangle motif and is
computed from the **pre-pruning** AC so it cannot be triggered by pruning itself. The two levers
are orthogonal in intent but **not** in code: a future edit to either must keep `exempt` computed
before `cap_order` is consumed, and must not assume `cap_order` is `range(N)`. There is a
detailed in-code comment block at `:1858-1888` saying exactly this; keep it in sync.

### Rejected alternatives

* **A canonical slot post-pass** (Lane 2's own going-in design) — refuted: this is not a
  slot-labelling problem, so no relabelling at that seam can reach it.
* **Fixing `_align_to_pai`'s `(i+1)**3` Z-moment** — the probe's leading suspect for `DUDREA`.
  Not the mechanism: `rotate` drift is 0, and the diagnosis lands in AC perception, upstream of
  alignment.
* **Relying on `OIN_CANONICAL_PERCEPTION`** — measured not to close this molecule (or the
  stereo class).
* **"Shortest contact wins"** as the capping rule — refuted by the 2.328 Å datum: the bond that
  flipped was third-shortest of six, so a distance rule would not have been deterministic here
  either.
* **Adding an atom-index tie-break to `_cap_key`** — would reintroduce the dependence.

---

## Dead ends, refutations, and costs accepted

### ⚠ THE REAL COST: on degenerate input, a LOUD FAILURE became a SILENT WRONG ANSWER

This is the cost that must not be lost, and it is the strongest argument the lever has against
it. `tests/fixtures/ticat3_generated_broken.xyz` is a deliberately broken fixture — a captured
TiCat3 *generated* structure whose clashing geometry over-connects a ligand fragment.

* **Old order:** the capping walk dead-ends, `get_lig_mol` cannot build the fragment, and
  `get_tmc_mol` raises a descriptive `ValueError` — **it fails loudly.**
* **New order:** capping highest-Z first lets the **titanium absorb the contested bonds**, so
  the walk does not dead-end and **perception SUCCEEDS — returning nonsense**: **48 atoms in 8
  fragments, seven bare `[H+]` ions, and a `[Ti-14]` centre.**

On a broken input, "perceives a nonsense graph" is worse than "fails loudly". Real data is
clean (`geometry_tag_shift` 0/298; 145 fixed with zero correctness regressions), so this is a
**degenerate-input concern, not a shipped-accuracy one** — but it is a real direction-of-failure
change and it is pinned rather than papered over.

**The follow-up worth considering is a sanity gate that rejects a perceived molecule containing
stranded bare-proton fragments. It was NOT added**, because charged hydrides are legitimate
species and such a gate needs its own corpus A/B before it can be trusted not to reject real
molecules.

**Consequence, and a lesson about contract tests.** `tests/unit/test_xyz2mol_errors.py` — the
unit pin for TASK-41/WS-1, that `get_tmc_mol` raises a descriptive `ValueError` instead of
returning a bare `None` on ligand-perception failure — **had to be rewritten as fault
injection**, because its fixture no longer triggers the error path. The contract still matters
(the original failure surfaced as the opaque *"cannot unpack non-iterable NoneType object"*,
because the sole convert-path caller in `core/translator.py` unpacks a 2-tuple), so it is now
exercised by making `get_lig_mol` fail directly:

```python
with mock.patch.object(xyz2mol_module, "get_lig_mol", return_value=(None, 0)):
    with self.assertRaises(ValueError) as ctx:
        get_tmc_mol(_FIXTURE, 0, with_stereo=False)
```

Note the injected value is `(None, 0)`, not a bare `None`: `get_lig_mol` returns a 2-tuple that
the call site unpacks *before* the `if not lig_mol` guard, so injecting a bare `None` would fail
in the unpack — which is the very `TypeError` this contract exists to prevent, raised from the
wrong place.

**The generalisable lesson, in the test's own words:** *"A contract test that depends on a
fixture staying unbuildable is one perception improvement away from silently testing nothing"* —
which is exactly what happened, and it would have hidden a regression back to the bare-`None`
`TypeError`.

A second test was added to pin the degenerate behaviour rather than pretend it is fine:
`test_broken_fixture_perceives_a_degenerate_graph_under_stable_metal_ac` asserts the fixture
yields more than one fragment and at least one bare-proton fragment. It is **not an
endorsement** — it is a tripwire: if a future sanity gate makes this raise again, or the capping
order changes, the test fails and the module docstring gets revisited instead of the behaviour
drifting unobserved.

### ⚠ A METHODOLOGY ERROR: comparing adjacency matrices without canonicalising atom order

An early A/B in this programme reported **"0/36 identical, 44–168 bond differences."** That
should have been an immediate tell: 44+ differing bonds on a 60-atom molecule would make it
unrecognisable. **The comparison was ORDER-SENSITIVE**, and the thing being compared writes its
atoms in a different order. Redone order-insensitively (bond element-pair multiset + (element,
degree) multiset):

| | result |
|---|---|
| atom **symbol sequence** identical | **0/36 (0%)** |
| **graph** fingerprint identical | **16/36 (44%)** |

**Comparing adjacency matrices without canonicalising atom order manufactures differences.** In
this specific instance the corrected reading also flipped the engineering conclusion twice over:
the graph *is* preserved 44% of the time, but the atom order **never** is — so the memo whose
hit-rate was being estimated cannot hit at all, because its key is a permutation away; and 56%
of generated structures do not preserve the input's graph even up to isomorphism, which is a
large number bearing on the `structural` bucket.

> **Provenance, stated so nobody hunts for it in the wrong place.** This particular 0/36 → 16/36
> correction is recorded in `docs/agentic-notes/v0.4.5/V045_STATUS_2026-07-25.md` §"I tested the lane's own follow-up
> hypothesis, and my first test was broken", and the test in question was the `encspeed` lane's
> question *"is a generated conformer's ligand AC byte-identical to the input's?"* — **not** an
> `OIN_STABLE_METAL_AC` A/B. It is recorded here because it is the same class of error as this
> lane's subject matter (atom-order sensitivity in AC comparison) and because any future A/B on
> this lever will be tempted to compare AC matrices directly. If you do, canonicalise first.

### Other refutations and accepted costs

* **The `(i+1)**3` Z-moment hypothesis** — plausible, written down, and wrong for this molecule.
  0 rotation drift across 225 × 2 trials is what killed it.
* **The residual is untouched by this lever.** After promotion, byte-stability is 69.6% and 16
  molecules still have an unstable comparison key. The residual `slot_renumber` class is
  structurally harder — **32/32** of Lane 2's residual pairs are `same_vcolor_identical`, so no
  relabelling at that seam can close them.
* **`slot_renumber` counts ROSE** (42 → 74) while byte-stability rose. That is
  **reclassification, not regression**: a molecule that previously drifted in *both* its ligand
  body and its slot numbers was counted under `rdkit_canonical` (first matching subclass), and
  once the body is canonical its remaining slot drift reclassifies. Lane 2 established **0
  molecules broken** in either arm with per-molecule accounting. Anyone reading the subclass
  histogram cold will misread this.
* **An earlier version of the promotion-gate document published 2 of 3 shards** as 62.0% → 73.5%
  and key-broken 39 → 8 (−79%). The third shard was harder, so absolute levels are lower and the
  key reduction is −73%, not −79%. The **delta held at exactly +11.5 points** across both reads.
  Partial absolutes should not stand as if final.
* **What the promotion does *not* buy:** it barely moves the round-trip headline, and that is
  structural — a canonicality defect lands in `key_equal`, which **already counts as passing**.
  Of the 332 closeable molecules in the gap to the ~98.45% ceiling, `hard_fail` is **306 (92%)**.
  What promotion buys is that the number becomes **meaningful**: before it, 13% of molecules
  encoded a different absolute stereochemistry depending on the input file's atom order and the
  key moved for 1 in 5, so "round-trip success" was partly a property of how the XYZ happened to
  be numbered.

---

## Where it landed

**Merged and released.** `git log --oneline main..swimlane/v045-lane2` is empty. Released on
local `main` as tag `v0.4.5` (`0d165845`); integration commit `1450b5ce` (`release(v0.4.5):
integrate 16 lanes and PROMOTE the six canonicality levers`). ⚠ Per standing project
instruction, v0.4.3/0.4.4/0.4.5 are **local-only and must not be pushed.**

| commit | subject |
|---|---|
| `8bf9df61` | `fix(v0.4.5): make AC valence-capping order-independent (OIN_STABLE_METAL_AC, OFF)` — the fix, 1 file, +41/−4 |
| `d493f192` | `tools(v0.4.5): geometry-tag shift check — the veto on OIN_STABLE_METAL_AC` |
| `552c06a6` | `docs(v0.4.5): OIN_STABLE_METAL_AC passes its geometry veto — recommend promotion` |
| `51354223` | `docs(v0.4.5): promotion gate PASSED — recommend all six canonicality levers ON` |
| `2a1c50df` | `docs(v0.4.5): replace the promotion-gate figures with the full-sample result` |
| `840eab84` | `lane2(v0.4.5): pin exactly which goldens the lever moves, and that each is key-identical` |
| `20044883` | `tools(v0.4.5): rotation/renumbering canonicality probe (replaces the failed trust gate)` |

**Final state**

| item | value |
|---|---|
| lever | `OIN_STABLE_METAL_AC` |
| default | **ON** — `src/oinsmiles/oin/levers.py::_DEFAULT_ON` |
| opt out | `OIN_STABLE_METAL_AC=0` (also `false`/`no`/`off`/empty) |
| code | `src/oinsmiles/utils/xyz2mol_local.py:1858-1913` (`xyz2AC_obabel`, `:1809`); comment block `:1858-1888` documents both defects in that loop |
| composes with | `OIN_BORON_CAGE` (`boron_cage_vertices`, `:1744`) in the same statement |
| instrument | `tools/geometry_tag_shift.py` (`--lever OIN_STABLE_METAL_AC --n 300 [--shard 1:4]`) |
| probe | `tools/canonicality_probe.py` (`--n 300 --trials 2`, seed 42 fixed so every arm samples the same molecules) |

**Guard tests**

* `tests/unit/test_xyz2mol_errors.py::TestXyz2MolErrors::test_get_tmc_mol_raises_valueerror_when_get_lig_mol_fails`
  — the error-path contract, now fault-injected.
* `tests/unit/test_xyz2mol_errors.py::TestXyz2MolErrors::test_broken_fixture_perceives_a_degenerate_graph_under_stable_metal_ac`
  — the degenerate-input tripwire (asserts >1 fragment and ≥1 bare-proton fragment).
* `tests/unit/test_regression_stability.py::TestRegressionStability` (6 golden fixtures:
  `test_cisplatin`, `test_transplatin`, `test_cis_ptcl2en`, `test_ferrocene`,
  `test_fac_irppy3`, `test_mer_irppy3`) — byte-identical on the opt-out path.
* `tests/unit/test_facmer_key.py` + `tests/unit/test_isomer_divergence.py` — the over-folding
  veto: fac ≠ mer and cis/trans stay distinct **raw and at key level** with all six levers ON.
* `tests/unit/test_stable_stereo_mirror.py` — the sibling lever's must-not-be-stable-because-constant
  guard, re-run under all six levers ON (10/10 mirrors differ).

**Reproduce**

```bash
cd /home/tjmustard/Documents/GitHub/oin-v045-trial
export PYTHONPATH=$PWD/src
V=/home/tjmustard/Documents/GitHub/OIN-SMILES/.venv/bin/python
DS=/home/tjmustard/Documents/GitHub/OIN-SMILES/tmCAT-tmPHOTO_xyz_dataset

# the geometry veto (the lever-specific evidence)
$V tools/geometry_tag_shift.py --lever OIN_STABLE_METAL_AC --n 300 --dataset "$DS" --out <dir>

# the three hand-picked worked examples
PYTHONPATH=src .venv/bin/python tools/canonicality_probe.py \
    --only FEQFIS_comp_0,DUDREA_comp_0,CEBVIR_comp_0 --trials 3 -v
```

---

## Open questions / for the next agent

1. **The stranded-bare-proton sanity gate.** The named, deliberately-unbuilt follow-up. A gate
   rejecting a perceived molecule that contains isolated bare-`[H+]` fragments would restore
   loud failure on degenerate inputs, but **charged hydrides are legitimate**, so it needs its
   own corpus A/B (how many currently-passing molecules contain a legitimate isolated hydride?)
   before it can ship. If you build it, `tests/unit/test_xyz2mol_errors.py`'s second test is the
   tripwire that will tell you it worked — and its module docstring must be updated with it.
2. **Reconcile the coordination-degree numbers** (5/4 in the narrative vs 6/7 in the
   measurement). See the flagged inconsistency above.
3. **The two residual renumbering classes.** `FEQFIS_comp_0`'s stereo flip and
   `CEBVIR_comp_0`'s aromaticity collapse survived all four v0.4.5 perception levers as
   measured; `OIN_STABLE_STEREO` later closed 8/10 of the stereo class, leaving 2 drifting (1 at
   key level). Lane 8's own note says trivalent phosphorus donors are unhandled, and the reason
   is **not** a Lane 8 defect: `stable_stereo.py:112-118` only *corrects tags that already
   exist*, so a P donor whose tag was cleared by the Zone-A rule (`core/chirality.py:722-727`)
   has nothing to restamp. Making metal-locked donor tags exist is P3's job.
4. **`CEBVIR`-style aromaticity collapse is still unowned.** It is consistent with `AC2BO`'s
   "arbitrary resonance form" and `get_UA_pairs`' non-unique `nx.max_weight_matching` being
   atom-order dependent — i.e. the *other* valence-order lane's territory (see
   `LANE-valence-search.md`), where the same matching non-uniqueness is measured on the over-cap
   population.
5. **`OIN_BORON_CAGE` is now also default-ON (v0.4.6).** The two levers' composition inside the
   capping loop was never A/B'd **together** on the corpus — each was measured against a
   lever-OFF baseline. A combined `geometry_tag_shift` run with both ON versus both OFF would
   close that gap cheaply.
6. **Is the invariant tie-break sufficient on highly symmetric molecules?** `_cap_key` leaves
   genuinely equivalent atoms tied, and `sorted` is stable — so among tied atoms the **input
   order still decides**. On a symmetric molecule where two equivalent atoms compete for the same
   bond, order dependence can therefore survive. Nothing measured this; the 0/8 result on DUDREA
   does not cover it. A targeted probe on a high-symmetry cage or a homoleptic complex would
   settle it.
7. **The comparison key's historical error term.** The key changed under renumbering for 21.1%
   of sampled molecules before this work and 5.4% after, so **every absolute accuracy figure
   this project reported before v0.4.5 carries an unaccounted systematic term.** Relative A/B
   results are unaffected (both arms share the input ordering). Do not re-derive old absolute
   numbers as if they were clean.
