# Injectivity Y3 residuals (Lane 7) — the torsion-aware oracle, the donor swap, the missing operators

**Status:** measurement wave — **no encoder output changed**. Parent:
`INJECTIVITY_Y3_UNKNOWN_UNKNOWNS.md`, whose Part 2 residue and Part 3 open questions this
closes. Reproducers: `tools/injectivity/{torsion_oracle,config_split,twin_operators,positional_isomers}.py`.

Wave 3 ended with three things unresolved and one thing unbuilt:

| left open | why it could not be answered then |
|---|---|
| 70 of 299 structures in an ambiguous residue | the oracle was a **rigid** superposition test, so it could not tell a flexible molecule's mirror-conformer from a genuine enantiomer |
| "does the key fold a donor swap?" — undetermined | the only probe was hand-written OIN strings, which never established the two strings named *different isomers* |
| `@SP` metal stereo | a mirror cannot produce a square-planar diastereomer, and no other operator existed |
| a Δ/Λ tris-chelate and a 4-different-donor square-planar fixture | never built; Lane 5 had only `fac-Ir(ppy)3` |

---

## Part 0 — the two fixtures Lane 5 is validated on

`fac-Ir(ppy)3` is the *easy* case for metal Δ/Λ: its three chelates are unsymmetric (C,N), so a
descriptor that actually encodes **fac/mer** rather than **helicity** would still appear to work
on it. That is the same shape of trap the Y2 wave hit, where a reflection-*invariant* axial
token passed every guard because the only fixture exercised the easy single-axis case.

Both fixtures are real crystal geometries, selected by a census of the 26 230-structure corpus
for their archetype (single metal, ≤ 90 atoms).

| fixture | what | oracle on its mirror | metal descriptor | vdW clashes |
|---|---|---|---|---|
| `tests/fixtures/ZUMNEC.xyz` | tris(catecholato)Mo, 37 atoms | **distinct, 1.34 Å** | `@OH` permutation **11 → 9** | 0 |
| `tests/fixtures/JEGKOW.xyz` | Rh(I) square planar, donors N/P/C(carbonyl)/I, 31 atoms | **not distinct, 0.099 Å** | `@SP` permutation 2 | 0 |

**ZUMNEC is the only homoleptic tris-bidentate in the corpus whose two donors per chelate are
symmetry-equivalent** (`CanonicalRankAtoms(breakTies=False)` gives them a single rank). That is
what makes it the guard: with equivalent donors there is no fac/mer distinction at all, so metal
helicity is its *sole* stereogenic element, and it carries no axial or bound-amine axis either.
Its 1.34 Å mirror RMSD is 2.7× the 0.5 Å threshold and 13–27× the achiral controls (0.05–0.10 Å).
It reproduces the P1 `key_blind` collapse today.

**JEGKOW's mirror is correctly *not* a new isomer, and demanding otherwise would be a category
error.** A square-planar complex is planar, so its coordination plane is a mirror plane and four
different donors give **diastereomers, not enantiomers**. Reflection is the wrong distinctness
operator for `@SP`; the right one is a donor swap — which is why Task C had to exist before this
fixture could be probed at all.

---

## Part 1 — the torsion-aware configurational oracle

### The instrument

`tools/injectivity/oracle.py::geometric_chirality` asks: *can the mirror be superimposed by a
proper rotation, over the graph automorphisms?* On a flexible molecule the answer is no even
when the two are the same isomer, because the mirror is a different **conformer**.

`tools/injectivity/torsion_oracle.py` asks the question that actually separates the two cases:

> **is the mirror reachable from this structure by rotating bonds?**

Rotating a dihedral is a continuous deformation; it cannot change any configuration — not a
tetrahedral centre, not the metal's coordination handedness. So the torsion orbit of the input
*is* the conformational orbit of its configuration. Reachable ⇒ conformational ⇒ **collapsing is
correct**. Unreachable ⇒ configurational.

**Why not a conformer pool.** Embedding the graph with `EmbedMultipleConfs` and asking whether
any conformer fits the mirror is unsound *here specifically*: the configurations under test
(metal Δ/Λ, atropisomerism) are exactly the ones RDKit does not carry in the graph, so a free
embed generates both handednesses and the pool fits the mirror every time. The test would answer
"conformational" unconditionally.

**The seed that makes the search work.** A reflection preserves every bond length and angle and
**negates every signed dihedral**. So the torsion vector that negates all rotatable dihedrals is
not a heuristic starting point, it is the physically correct guess, and a purely conformational
chirality is found there immediately. Random restarts alone are not a substitute: reaching the
mirror generally needs several torsions to flip *together*, and one-at-a-time coordinate descent
must climb before it descends.

**"No match" is evidence, not proof — so every verdict is paired with a control.** Same molecule,
same optimiser, same budget, recovering a randomly-torsioned copy of the structure that is
reachable *by construction*.

| control | mirror | verdict |
|---|---|---|
| converges | converges | `conformational` |
| converges | does not | `configurational` |
| does not converge | — | `inconclusive` — no evidence either way |

`TorsionVerdict` carries `budget`, `threshold`, `d_control` and the automorphism counts, so a
conclusion cannot be quoted without the search that produced it.

### Cross-validation — ten known answers, all correct

| structure | expected | `d_mirror` | `d_control` | verdict |
|---|---|---:|---:|---|
| `EDOQIZ` — unsubstituted biphenyl on linear Au | conformational | 0.076 | 0.075 | ✅ conformational |
| `WAVGOS` — bis-NHC, propargyl arms | conformational | 0.048 | 0.025 | ✅ conformational |
| `PERPIO` — Ni alkynyl + phosphine | conformational | 0.292 | 0.115 | ✅ conformational |
| `PdCl2-R-BINAP` — hindered biaryl | configurational | 3.448 | 0.032 | ✅ configurational |
| `YESKOZ` — two hindered axes | configurational | 2.819 | 0.017 | ✅ configurational |
| `fac-Ir(ppy)3` — metal Δ/Λ | configurational | 2.673 | 0.000 | ✅ configurational (rigid) |
| `ZUMNEC` — metal Δ/Λ | configurational | 1.065 | 0.000 | ✅ configurational (rigid) |
| `POJJOP` — metal-bound amine | configurational | 1.707 | 0.009 | ✅ configurational |
| `CisPlatin` | achiral | 0.003 | — | ✅ rigid achiral |
| `JEGKOW` | achiral | 0.067 | — | ✅ rigid achiral |

`WAVGOS`, `PERPIO` and `EDOQIZ` are the three cases Wave 3 triaged **by hand** and called
conformational. The tool reproduces all three, which is the closest thing available to an
external check.

### Two bugs the cross-validation caught

Both made the tool confidently wrong rather than merely imprecise, which is why the
cross-validation set had to include cases in *both* directions.

**1. Cutting the atom instead of the bond.** Rotatability was tested by asking whether the graph
disconnects when an *atom* is blocked. A metal is a cut vertex, so blocking it detaches every
ligand from every other, and each metal–donor bond of a **chelate** looked rotatable. The result:
`fac-Ir(ppy)3`'s Δ/Λ mirror was "reached" at `d_mirror = 0.042 Å` by swinging whole ligands off
their own chelate rings. A confident wrong answer, not a near miss.

**2. Automorphism starvation.** Automorphisms were enumerated on the H-explicit graph, where
methyl rotations consume the whole `maxMatches` budget on permutations that leave every heavy
atom fixed. On `EDOQIZ` (two *tert*-butyls) 4 000 full-graph matches collapse to **6** distinct
heavy images; the heavy skeleton has **864**. Since more automorphisms can only *lower* an RMSD,
starvation inflates every number — and that is why a freely rotating biphenyl read as
configurational.

> ### ⚠ `oracle.py` has the same starvation, and it is left in place
>
> The rigid oracle enumerates on the H-explicit graph too. It changes **no curated fixture's
> verdict** (fac-Ir(ppy)3, BINAP, POJJOP and CisPlatin all stay on the same side of the
> threshold), so nothing published in `BASELINE.md` §3 moves. But it inflates the rigid mirror
> RMSD on methyl-rich species — `EDOQIZ` reads 2.55 Å starved and 0.50 Å complete — which means
> **the rigid oracle over-reports chirality at dataset scale**, and the `uu_hunt` "chiral" gate
> inherits that. Correcting it would change the published per-fixture RMSDs, so it is flagged
> here for the release owner rather than changed inside a measurement lane. The residue split
> below is unaffected: `config_split` uses the complete enumeration, so it re-decides every case.

### ⚠ What `conformational` does *not* mean

The search models **geometry, not energy**. A rotational barrier is invisible to it, so an
isolable atropisomer whose axis is a free (non-ring-locked) torsion is "reachable by rotating
bonds" and reads `conformational`. That is the strict geometric truth and the chemically
incomplete one.

This does not corrupt the residue split, because the barrier question is triaged *upstream*:
`uu_hunt` routes any structure with a stereogenic hindered axis into `known_axis` (P2) before the
residue is formed, using `oin.axial.detect_axial_axes`' hindrance heuristic. The two instruments
are complementary — `detect_axial_axes` supplies the barrier judgement, `torsion_oracle` supplies
the geometry — and neither should be used for the other's question.

---

## Part 2 — the donor swap

See `results-injectivity-y3/positional_isomers.md` for the generated tables.

### Why the original probe proved nothing

Hand-writing two OIN strings with donors in swapped slots makes the key look blind. It is not
evidence: for a square-planar complex with two identical ancillary ligands, swapping slots 0↔1 is
a reflection of a **planar, hence achiral** complex, so the two strings name the same molecule
and folding them is *correct*. The probe never established its two strings denoted different
isomers.

### Two independent geometry-driven lines

**Line 1 — real corpus pairs.** Two crystal structures with the same constitution (identical
canonical SMILES of the perceived complex) and a different **trans-donor multiset**. For a fixed
constitution that multiset is a configurational invariant — unchanged by any rotation of the
complex or of any bond — so two members whose multisets differ cannot be the same isomer. No OIN
string enters that certification.

The first "pair" the scan produced was a **false positive worth recording**: `BOCYEA_comp_0` vs
`BOCYEA_comp_1`, two crystallographically independent copies of the *same* compound, separated
only by one donor–M–donor angle sitting either side of the 150° trans cutoff. The scan now
requires the arrangement to differ at **every** cutoff in 140/150/160°, requires the trans-pair
*count* to match (a differing count means a straddled boundary or a different coordination
number, neither of which is a positional isomerism), and excludes same-refcode pairs.

**Line 2 — `swap_donor` on real geometry.** Symmetry-equivalent donor pairs are skipped —
exchanging those is the identity on the isomer, which is precisely the mistake the hand-written
probe made. Swaps of two **trans** donors are the built-in negative control: a trans exchange is a
180° rotation of the whole complex, so it is the same isomer and must read `distinct=False`.

---

## Part 3 — the missing twin operators

`twin_collision` shipped one operator, `mirror_z`. It is trivially valid (a coordinate transform
cannot break a structure) and answers exactly one question: whole-molecule enantiomerism.

| operator | what it reaches that a mirror cannot |
|---|---|
| `swap_donor` | **diastereomerism** — a square-planar complex is achiral, so only a donor exchange distinguishes `@SP1`/`@SP2`/`@SP3` |
| `invert_axial` | **one axis at a time** — a mirror flips every axis together and cannot localise the loss |
| `invert_tetrahedral` | a single tetrahedral centre in place |

These are **structural edits**, so they can produce nonsense. Every twin carries its vdW clash
count before and after, `probe_operator` refuses to score a twin the gate rejects, and
`invert_tetrahedral` generates every candidate exchange, ranks by clash and relaxes the winner in
torsion space before the gate rules.

### `invert_axial` — the Lane 4 result

`YESKOZ` carries two symmetry-equivalent hindered axes of opposite sign. Flipping **one** turns
the (+,−) diastereomer into (−,−) — an edit no mirror can make.

| build | oracle | raw strings | key | verdict |
|---|---|---|---|---|
| default (`OIN_EMIT_AXIAL` off) | distinct | **identical** | equal | `encoder_blind` |
| `OIN_EMIT_AXIAL=1` | distinct | **differ** | equal | `key_blind` |

**The Y2 axial token works on a multi-axis molecule.** The Wave-2 cohort A/B could not establish
this because its multi-axis arm failed for a *generator* reason (the torsions relaxed out of the
hindered window). This probe is generator-free, so that confound does not apply: with the token
on, the encoder separates the two diastereomers, and the residual fold now lives in
`compare.py`'s `_AXIAL_TOKEN_RE`, not in the encoder. → Lane 4.

### `invert_tetrahedral` — measured scope limit

The inversion is correct: POJJOP's signed tetrahedral volume flips sign. But **every stereocentre
in the fixture set is locked** — inside a chelate ring (BDPP, DPDME) or bound to the metal
(POJJOP) — and for all of them the rigid exchange drives a substituent into the coordination
sphere. Torsion relaxation reduces the damage (POJJOP 3 → 2 clashes, DPDME 10 → 4) but does not
clear it. Reaching those would need bond-angle relaxation, which a generator-free instrument must
not do. For a locked centre, `mirror_z` or `invert_axial` remains the instrument; this operator
is for stereocentres on a freely rotating pendant.

---

## Reproduce

```bash
PYTHONPATH=$PWD/src python -m tools.injectivity.report --probes --operators
PYTHONPATH=$PWD/src python -m tools.injectivity.torsion_oracle tests/fixtures/PdCl2-R-BINAP.xyz
PYTHONPATH=$PWD/src python -m tools.injectivity.uu_hunt --n 300      # writes the residue
PYTHONPATH=$PWD/src python -m tools.injectivity.config_split         # splits it
PYTHONPATH=$PWD/src python -m tools.injectivity.positional_isomers --scan
```

Output lands in `results-injectivity-y3/` (gitignored, regenerable, seeded at 42).
