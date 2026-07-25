# Lane 3 — the eta winding residual, characterized

Status: v0.4.5. Instruments: `tools/winding_diag.py`, `tools/eta_core_chirality.py`,
`tools/reencode_ab.py` + `tools/roundtrip_bucket_report.py`, `tools/canonicality_probe.py`.

> **Headline.** The lane's stated premise — "the 13 residuals come from the geometric
> fallback tiers" — is **wrong for all but one of them**. The encoder's winding sign is
> already invariant under rigid rotation, and in 5 of the 6 cases that survive on current
> `main` the geometric heading tier never fires at all. One case is a genuine
> generator-produced **enantiomer** that must NOT be folded.

---

## 1. The published list of 13 is stale

`winding_star_drift = 13` comes from the **v0.4.2** capstone sweep
(`results-capstone-v042/bucket_report.json`). Re-encoding those same 13 stored geometry
pairs with current `main` reclassifies them:

| current bucket | n | molecules |
|---|---:|---|
| `key_equal` / **`winding_star_drift`** | **6** | GAMJAG, GIPDEQ, QIGZAJ, SIMDIE, TEYXEA, WUHRIB |
| `key_equal` / `rdkit_canonical` | 3 | BOJMAQ, MOHKOL, SATKIJ |
| `structural` | 3 | ABETIK, IFICAD, OVUBEO |
| `facmer_divergent` | 1 | FAHCIC |

v0.4.3/v0.4.4 already moved 7 of the 13 out of the class. **The lane's true worklist is 6.**
(BOJMAQ is not a winding case at all — its two eta groups keep identical per-slot windings
and merely trade slot numbers, i.e. Lane 2. MOHKOL/SATKIJ do show the sign pattern below,
but with body drift on top, so Lane 1 owns them first.)

## 2. Per-molecule table

`tier` = the tier that chose the heading atom, from `tools/winding_diag.py`, measured
separately on the input and on the stored generated structure. `rot-inv` = the emitted
string is identical across N random **proper** rotations of the input.

| molecule | symptom | deciding tier (input / generated) | rot-inv | `generated == mirror(input)` | root cause | verdict |
|---|---|---|:---:|:---:|---|---|
| **GIPDEQ** | star moves (B → CH₂), sign unchanged | **2-GEOMETRIC / 2-GEOMETRIC** | yes | no | `_canonical_heading_atom` returns `None` (boron ylide `[CH2]=c1cc(C)cc(C)c1=B(...)` will not kekulize), so Tier 1 fell through to the geometric heading, which tracks the embedding | **fixed** (Fix A) |
| **GAMJAG** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | two **separate, identical** benzylindenyl fragments at interchangeable TET slots; the slot-labelling maximization is a genuine 2-way tie broken by the Kabsch fit | **fixed** (Fix C) |
| **SIMDIE** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | symmetric Me₂Si-bridged bis(2-Me-4-Ar-7-OMe-indenyl): the **meso** diastereomer, whose two mirror-related spellings are one achiral compound | **fixed** (Fix C) |
| **TEYXEA** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | as SIMDIE (Si-bridged, two equivalent rings) | **fixed** (Fix C) |
| **WUHRIB** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | as SIMDIE (butenyl-bridged bis(tBu-Cp), two equivalent rings) | **fixed** (Fix C) |
| **QIGZAJ** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | Cp ↔ fluorenyl are **not** automorphic, so the two spellings are **genuine enantiomers**. The generator built the wrong one | **not a bug — must not be folded** |

### Three distinct causes, not one

1. **Geometric heading fallback** (GIPDEQ, 1 case) — the lane's stated premise. Real, and
   the only case where a geometric tier is load-bearing.
2. **Eta slot-labelling / automorphism freedom** (GAMJAG + SIMDIE/TEYXEA/WUHRIB, 4 cases) —
   the encoder is correct and rotation-invariant; it simply has two equally valid spellings
   for one compound and picks between them by geometry.
3. **A generator stereochemistry error the key over-folds** (QIGZAJ, 1 case) — see §5.

## 3. What was measured, and what it rules out

- **The sign math is already rotation-invariant.** `signed_circulation` computes
  `cross(v_star, v_next) · axis` with `axis` the ring's *own* metal→centroid vector
  (`_determine_winding`, `oin_aligner.py`). All three factors rotate together, so the sign
  is invariant by construction. Confirmed empirically: N random proper rotations of every
  case give **one** distinct string. The orchestrator's independent
  `canonicality_probe.py` smoke found the same thing corpus-wide (zero `rotate` drift).
- **Therefore `winding_star_drift` is invisible to `canonicality_probe.py`.** That probe
  holds the structure fixed and varies presentation. These residuals are *input structure
  vs a differently-generated structure of the same compound* — a different question. The
  acceptance number for this lane must come from the `reencode_ab` + bucket-report
  `key_equal` sub-split, which the orchestrator confirms is the trustworthy part of that
  instrument.

## 4. Why the mirror test is the right question, and why RMSD cannot answer it

Every one of the 6 satisfies `generated == mirror(input)` **byte for byte**. So the question
is not "did the encoder drift?" but "**is `mirror(M)` the same compound as `M`?**"

`tools/injectivity/oracle.py::is_distinct_enantiomer` answers "chiral" for all 6 (RMSD
2.6–4.4 Å, automorphism cap hit) — but that is **conformational** chirality: these molecules
have pendant aryls, tBu rotors and, for the unbridged ones, freely rotating rings, so a
genuinely achiral meso complex still fails to superimpose on its mirror.
`tools/eta_core_chirality.py` was built to strip that away (metal + coordinated ring systems
+ bridge + one substituent shell) and **does** calibrate correctly on rigid ferrocene
(ACHIRAL, 0.103 Å, 200 core automorphisms — after flattening formal charge, which
perception localizes onto one arbitrary Cp carbon and which otherwise collapses the
automorphism group to 8 and reports ferrocene as chiral at 0.81 Å). But it still cannot
decide a metallocene whose rings rotate freely about the metal axis, so it is **not** the
oracle for this lane. It is kept, with that limitation documented, because it is a correct
instrument for rigid cores.

### The criterion that does work — and is exact

Two OIN spellings denote the same compound exactly when related by a **proper rotation of
the coordination polyhedron** composed with a **graph automorphism of the ligand assembly**.
Both operations carry each ring's winding character along *unchanged*: a proper rotation
preserves handedness, and an automorphism mapping one eta ring onto another preserves the
cyclic sense (an automorphism that *reverses* it is exactly the orientation-free case,
already forced to a fixed `>`). So `(ring identity, winding char)` travels as a unit, and
the only freedom is **which interchangeable slot each unit is written at**.

This reproduces the textbook chemistry as a consistency check:

| arrangement | signs | mirror | verdict |
|---|---|---|---|
| **rac** ansa-metallocene | `(>,>)` | `(<,<)` — different multiset | chiral, two distinct strings ✅ |
| **meso** ansa-metallocene | `(>,<)` | `(<,>)` — same multiset, reachable by the C2 that swaps the two eta slots | achiral, must be one string ✅ |
| Cp/fluorenyl (QIGZAJ) | `(>,<)` | `(<,>)` — rings **not** automorphic, so not reachable | chiral, two distinct strings ✅ |

## 5. QIGZAJ: the key over-folds an enantiomer pair (a finding for Lane 2)

QIGZAJ buckets as `key_equal` — the key says input and generated are the same isomer. They
are not. `compare.py::_parse_vertex_colors` colours every slot marker in a fragment with
that **whole fragment's** body, so for an ansa ligand the Cp slot and the fluorenyl slot
receive the *same* colour. The tetrahedral C2 that swaps slots 0↔1 then looks
colour-preserving and `_polyhedron_signature` folds the two spellings. Per-eta-ring identity
(as in `_eta_automorphism_class`) would fix it. **Lane 3 deliberately does not fold QIGZAJ**,
so it remains a visible `winding_star_drift` — correctly.

## 6. The fixes

### Fix A — eliminate the geometric heading tier (default ON, `oin_aligner.py` step 4a)

Both branches of the content-canonical tier now call `_topological_heading_atom`, which has
its own three sub-tiers (strict canonical rank → valence-tolerant symmetry graph → lowest
constituent index), *all* functions of the graph alone. Previously the non-orientation-free
branch called `_canonical_heading_atom` directly and, on `None`, fell through to the
geometric loop. The change is scoped to fragments where the strict rank is unavailable —
a population that is by construction already embedding-unstable — so it cannot destabilize
a pair that is byte-exact today: identical strings imply identical fragment graphs, and the
topological heading is a function of the graph.

This follows the lane's instruction to *eliminate* a geometric fallback rather than make it
deterministic. After it, the geometric heading loop is dead for eta groups; it is retained
as a fail-safe.

### Fix C — canonical winding across automorphic eta groups (`OIN_CANONICAL_ETA_WINDING`, default **OFF**)

`OINDiscreteAligner._canonical_eta_winding` implements §4's criterion:

1. collect the eta groups whose winding is load-bearing (skipping orientation-free rings);
2. colour every occupied slot by `(chem_id, eta automorphism class)` — **winding excluded**,
   because we are deciding which slots may be *relabelled* and that must not depend on the
   characters about to be reassigned;
3. keep the proper rotations (`_brute_force_symmetries`, all from `Rotation.from_euler` —
   no reflection can enter) that preserve every slot colour;
4. within each orbit of eta slots that are **all in one automorphism class**, sort the
   winding characters and reassign in ascending slot order;
5. **guard:** only when the induced group realizes *every* rearrangement of the orbit
   (`|induced| == |orbit|!`). A rotation group may offer only cyclic permutations on a
   3-orbit, and assigning a sorted sequence the group cannot reach would fold genuinely
   distinct isomers. Otherwise fail safe to today's behaviour.

`_eta_automorphism_class` is the canonical SMILES of the fragment with every constituent
atom given the same atom-map number — two eta groups match exactly when some automorphism
of the fragment carries one constituent set onto the other. It works across fragments (two
identical separate indenyls) and within one (a bridged pair), and it deliberately skips
sanitization so an unkekulizable borate still gets an id.

**Why this cannot destroy stereochemistry.** No reflection is applied. On a rac
diastereomer both rings carry the same character, so the sort is a no-op and the mirror
stays a different string. Measured on 8 same-sign corpus cases (MECJOU, SERTUE, WOHNAJ,
EGIBEB, ODUFUO, HOHGEQ, RERQIO, XIBQEE): the lever changes **nothing**, and
`mirror != input` holds with it both OFF and ON. This is the guard the v0.4.4 axial wave
lacked when it sorted a token by sign and silently made it reflection-invariant.

## 7. Residual

`winding_star_drift` goes 6 → **1**, and the 1 is QIGZAJ, which must not be fixed here.
Closing it means teaching the *key* per-eta-ring colour (Lane 2) and fixing the
**generator**, which built the wrong enantiomer — not the encoder.
