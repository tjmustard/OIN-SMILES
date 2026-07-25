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
| **GAMJAG** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | two **separate, identical** benzylindenyl fragments; ε = +1, signs opposite ⇒ **achiral**, so the two spellings are one compound | **fixed** (Fix C) |
| **SIMDIE** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | Me₂Si-bridged bis(2-Me-4-Ar-7-OMe-indenyl); ε = −1, signs opposite ⇒ **chiral** — input and generated are **enantiomers** | **not an encoder bug** |
| **TEYXEA** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | as SIMDIE (ε = −1, opposite ⇒ chiral) | **not an encoder bug** |
| **WUHRIB** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | as SIMDIE (ε = −1, opposite ⇒ chiral) | **not an encoder bug** |
| **QIGZAJ** | both signs swap | 1a-strict / 1a-strict | yes | **yes** | Cp ↔ fluorenyl are **not** automorphic (ε undefined) ⇒ **genuine enantiomers** | **not an encoder bug** |

Plus one **latent** case the characterization surfaced that was not in the drift list at
all: **MECJOU** (ε = −1, signs *same* ⇒ achiral) encodes differently from its own mirror
image. Input and generated happened to agree, so no sweep ever flagged it — but it is the
same defect as GAMJAG and Fix C closes it too.

### Three distinct causes, not one

1. **Geometric heading fallback** (GIPDEQ, 1 case) — the lane's stated premise. Real, and
   the only case where a geometric tier is load-bearing.
2. **Un-canonicalized reflection freedom on an ACHIRAL eta arrangement** (GAMJAG, plus
   latent MECJOU) — the encoder is correct and rotation-invariant; it simply has two
   equally valid spellings for one compound and picks by geometry.
3. **A generator stereochemistry error that the key over-folds** (SIMDIE, TEYXEA, WUHRIB,
   QIGZAJ — 4 of 6) — the encoder is *right*; the generated structure is the enantiomer.
   See §5.

**So the majority of this lane's residual is not an encoder defect at all.**

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

### ⚠ The sense factor ε — where a first cut of this got it exactly backwards

The winding character is read off **each ring's own canonical SMILES order**. So an
automorphism carrying ring A onto ring B does *not* simply hand the character over — it
hands it over **possibly reversed**, depending on whether it maps A's cyclic order onto B's
in the same rotational sense (`ε = +1`) or the opposite one (`ε = −1`).

Applying the slot-swapping proper rotation together with the automorphism maps `(w₀, w₁)`
to `(ε·w₁, ε·w₀)`. That composite equals the mirror spelling `(−w₀, −w₁)` exactly when:

| ε | achiral arrangement is | typical structure |
|---:|---|---|
| **+1** | **opposite** signs `(>,<)` | two separate copies of one fragment — canonically ordered identically |
| **−1** | **same** signs `(>,>)` | two rings inside ONE bridged fragment — the canonical SMILES runs out along one ring and back along the other |

**ε is not a constant.** Measured: `+1` for the TiCat3/TiCat4 ligand, `−1` for SIMDIE and
MECJOU. A first implementation assumed the mapping and simply *sorted* the two characters;
that made the encoder **reflection-invariant** for every bridged case and folded a rac pair
— the v0.4.4 axial failure, reproduced exactly. It passed every guard written against the
easy fixture. What caught it was building the independent oracle and finding it reported
ACHIRAL for the same-sign case and CHIRAL for the opposite-sign one, the inverse of the
assumption. `_eta_swap_sense` now **computes** ε from the real fragment automorphism.

Verdicts from the computed ε agree with the independent geometric oracle **6/6**:

| case | signs | ε | ε-verdict | oracle |
|---|---|---:|---|---|
| TiCat3 | `(>,>)` | +1 | CHIRAL | CHIRAL (2.14 Å) |
| TiCat4 | `(<,>)` | +1 | ACHIRAL | ACHIRAL (0.03 Å) |
| SIMDIE | `(<,>)` | −1 | CHIRAL | CHIRAL (2.78 Å) |
| MECJOU | `(>,>)` | −1 | ACHIRAL | ACHIRAL (0.04 Å) |
| TEYXEA | `(<,>)` | −1 | CHIRAL | CHIRAL (2.65 Å) |
| WUHRIB | `(>,<)` | −1 | CHIRAL | CHIRAL (1.45 Å) |
| QIGZAJ | `(>,<)` | — (not automorphic) | no fold | CHIRAL (2.70 Å) |

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

`OINDiscreteAligner._canonical_eta_winding`:

1. collect the eta groups whose winding is load-bearing (skipping orientation-free rings);
2. colour every occupied slot by `(chem_id, eta automorphism class)` — **winding excluded**,
   because we are deciding which slots may be *relabelled* and that must not depend on the
   characters about to be reassigned;
3. keep the proper rotations (`_brute_force_symmetries`, all from `Rotation.from_euler` —
   no reflection can enter) that preserve every slot colour;
4. for each **2-orbit** of eta slots in one automorphism class, compute ε with
   `_eta_swap_sense` and apply the achirality test above. If achiral, emit the
   lexicographically smaller of the spelling and its mirror; if chiral, **leave it alone**;
5. **fail safe** on anything else — a larger orbit, a missing automorphism, an
   uncomputable ε — to today's behaviour.

`_eta_automorphism_class` is the canonical SMILES of the fragment with every constituent
atom given the same atom-map number — two eta groups match exactly when some automorphism
of the fragment carries one constituent set onto the other. It works across fragments (two
identical separate indenyls) and within one (a bridged pair), and it deliberately skips
sanitization so an unkekulizable borate still gets an id.

**Why this cannot destroy stereochemistry.** Only proper rotations and graph automorphisms
are applied, and the fold is gated on a *measured* achirality test rather than an assumed
one. Guards in `tests/unit/test_winding_canonical.py`: TiCat3/TiCat4 — a real rac/meso pair
of one Me₂Si-bridged bis(indenyl) ligand, differing *only* in eta winding — stay distinct
with the lever on; a chiral eta structure still differs from its mirror; an achiral one
folds only with the lever on; orientation-free metallocenes are untouched.

## 7. Residual — and why it is not zero

`winding_star_drift` goes **6 → 2** (GIPDEQ by Fix A, GAMJAG by Fix C; latent MECJOU also
closed). The remaining SIMDIE, TEYXEA, WUHRIB and QIGZAJ are **not encoder defects** — in
each the generator produced the enantiomer and the encoder correctly said so. Driving this
class to 0 in the encoder would mean folding enantiomers.

Closing them properly needs work in two other places:

- **the generator** — it does not reproduce the coordinated face of a bridged eta ring, so
  it returns the wrong diastereomer/enantiomer (the same shape of failure Lane 4 found for
  multi-axis atropisomers);
- **the key** (Lane 2) — it should stop folding them. `compare.py::_parse_vertex_colors`
  colours every slot in a fragment with the *whole fragment's* body, so an ansa ligand's
  two eta slots get the same colour and the tetrahedral C2 looks colour-preserving.
  Per-eta-ring colour, as in `_eta_automorphism_class`, would separate them.
