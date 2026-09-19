# Census C1 — the generated structure, judged without the encoder

**Date:** 2026-09-17 · **Branch:** `research/census` · **Plan:** [`CENSUS_PLAN.md`](CENSUS_PLAN.md)
**Instrument:** `tools/census/neutral_graph.py` (the ruler), `tools/census/g_vs_input.py` (runner +
controls), `tools/census/mirror_probe.py` (encoder-only follow-up).
**Population:** the 4,748 of 5,000 cohort molecules with a `*_generated.xyz` in
`results-v0.4.14-sweep/structures/`. Offline; the whole run is ~5 s on 10 cores.

## Headline

1. **`slot_renumber` (252 mol / 5.04 pts) is an ENCODER fault, not a generator fault.** 234 of
   the 252 generated structures are graph-identical to their input **and the molecule has no
   configurational handedness at all**. The roadmap files 201 of these as "G built the
   enantiomer" and chartered v0.4.17 L1 (CONSTRUCTION, 4.02 pts) on it. There is no enantiomer
   to build: `E` gives one achiral molecule two strings, differing only by a swap of slot labels
   between symmetry-equivalent donors of a symmetric ligand.
2. **68% of `structural` is confirmed generator construction by an encoder-free instrument**
   (331/484 graph ≠), and **29% is not** (142/484: graph identical, stereo identical or absent —
   the generator was right and the re-encode read it differently).
3. **Up to 380 of the 3,858 passes (9.8%) are flagged**; 40 are firm (mirror-image metal sphere,
   byte-identical string — the known P1 Δ/Λ blind spot, now counted at corpus scale).
4. **Side finding, invisible to the round trip:** UFF-tier output has metal–ligand bonds **15%
   too short** (median 0.818 × Σr_cov vs 0.965 in the crystal structures; Pt–Cl 1.90 Å vs 2.26–2.36).

## The ruler, in one paragraph

Distances only. No perception, no bond orders, no charges, no `oinsmiles` import (so it cannot
inherit the encoder's blind spots or run against the wrong source tree). Heavy-atom graph from
covalent-radius cutoffs; each metal sphere read three ways (firm / standard / generous) because
one fixed cutoff is the wrong model for a semi-coordinated or chelate-backbone atom; metal
spheres compared by proper vs improper superposition over every donor matching that preserves
symmetry class and chelate linkage; tetrahedral-type centres by signed volume, compared per
symmetry class plus a reflection-invariant relative-configuration check. Everything goes
through symmetry classes, never raw indices. Stated blind spots: E/Z, atropisomers, planar
chirality of a substituted η ring, ring cis/trans with symmetric paths — each can only cause a
false `SAME`, never a false `MIRROR`.

## Controls (all n = 5,000; every control structure is randomly rotated AND atom-permuted)

| Control | Construction | Must read | Result |
|---|---|---|---|
| identity | the input itself | ISO + SAME/NO_STEREO | **0 failures** |
| mirror | z-reflected input | MIRROR if ruler-chiral, else SAME | **0 failures** (1,771 = 35.4% ruler-chiral) |
| delete | smallest bound ligand removed | not ISO | **0 failures** |
| detach | largest bound ligand moved 3 Å out | not ISO | 24 false ISO (**0.48%**) |
| squeeze | every ligand slid 0.35 Å toward the metal | still ISO | 58 false non-ISO (**1.16%**) |

`identity`/`mirror`/`delete` are exact. `detach` and `squeeze` are **measured error rates, not
passes**: all 24 detach misses used the *firm* reading of the input, i.e. the ligand was already
only marginally attached in the crystal structure; the squeeze misses are inter-ligand clashes
the rigid slide itself creates. Read every count below with a ~1% floor in both directions.

What the controls do **not** establish: the ruler's *coverage* of chirality. `mirror` proves the
elements it sees negate correctly under reflection and are permutation-invariant; it cannot
prove no chiral molecule is read as `NO_STEREO`. That is checked from outside, below.

Three versions of the ruler were wrong before this one, and each was caught by a cross-check
rather than by its own controls — recorded because it is the project's recurring lesson:
(a) a single metal cutoff called 31% of passes non-isomorphic (cross-check: stored
`coordination` verdicts → 845/930 were *gained* second-shell contacts); (b) `squeeze` did not
exist until (a) showed the controls had no tolerance test; (c) the first sphere-chirality test
(unique-donor triples + skew lines) was incomplete (cross-check: the roadmap's 201
"enantiomers" read `NO_STEREO`, forcing the question of which instrument was blind).

## Result: ruler verdict × honest round-trip bucket

`ok:*` = graph isomorphic **and** stereo same (`same`) or no stereo element (`achiral`).

| bucket | n | ok:same | ok:achiral | MIRROR | DIASTEREO | G:sphere | G:ligand | no structure |
|---|---|---|---|---|---|---|---|---|
| byte_exact | 3858 | 1099 | 2379 | 40 | 104 | 159 | 77 | 0 |
| structural | 484 | 36 | 106 | 4 | 7 | 318 | 13 | 0 |
| hard_fail | 266 | 2 | 10 | 0 | 2 | 8 | 4 | 240 |
| key_equal / slot_renumber | 252 | 2 | **234** | 3 | 4 | 7 | 2 | 0 |
| key_equal / rdkit_canonical | 113 | 15 | 32 | 4 | 0 | 58 | 4 | 0 |
| facmer_divergent | 15 | 0 | 2 | 0 | 0 | 12 | 1 | 0 |
| encode_fail | 12 | 0 | 0 | 0 | 0 | 0 | 0 | 12 |
| **total** | 5000 | 1154 | 2763 | 51 | 117 | 562 | 101 | 252 |

Graph quality of the 4,085 isomorphic verdicts: 3,314 clean, 281 marginal (matched only under a
non-standard sphere reading), 490 `ISO_CLASH` (matched only with the covalent cutoff tightened to
0.88 — a non-bonded heavy-atom pair sits inside bonding range in the generated geometry).

Cross-checks against existing instruments: `structural` G:sphere **318** vs `attach_class_audit`
DETACHED **301** ✓ agree. `slot_renumber` MIRROR **3** vs `veto_residue_chirality`
MIRROR_MATCH **201** ✗ disagree — resolved in the next section.

## The disagreement, resolved: `E(mirror(x)) ≠ E(x)` for achiral `x`

`veto_residue_chirality.py` measured that the re-encoded string equals the encoding of the
input's mirror image. That measurement stands. The inference "therefore G built the
enantiomer" needs one more fact — *is the input chiral?* — and the ruler says no for 234 of 252.
By eye, 7 of 7 sampled agree, and the strings show the mechanism directly:

```
BIXTOU  in : …c1ccc(-c2cccc(-c3ccccn{4}3)n{3}2)n{5}c1     Ru(terpy)(N^N)Cl, Cs
        out: …c1ccc(-c2cccc(-c3ccccn{5}3)n{3}2)n{4}c1     outer terpy N: 4 <-> 5
BELYUP  phen N{3} <-> N{5}, fac-Re(CO)3      AGAVIQ  POCOP pincer P{3} <-> P{4}
AXESAY  DPEphos P{2} <-> P{3}               BOFLIT  salen N{0}<->N{2}, O{1}<->O{3}
```

The relabeling is induced by a graph automorphism of the whole complex, so both strings name the
same molecule. Encoder-only confirmation (`mirror_probe.py`, no generator involved):

| Group | n | `E(mirror) ≠ E(x)` |
|---|---|---|
| ruler-achiral, `slot_renumber` FAIL | 244 | **239 (98.0%)** |
| ruler-achiral, `byte_exact` PASS (random 600) | 600 | 58 (9.7%) |
| ruler-CHIRAL, `byte_exact` PASS (control — must differ) | 200 | 179 (89.5%) |

So ≈ 230 passing + 239 failing ≈ **470 molecules (9.4% of the corpus) are achiral yet
reflection-sensitive in `E`, and 51% of them fail** — a coin flip on which labelling the
round trip happens to land. (The coin flip alone does not prove achirality; a truly chiral
population with a hand-blind generator would show it too. The achirality rests on the ruler +
inspection.) This is also why `OIN_ACCEPT_STRING_EXACT` recovered 0.5% here against 28.7%
elsewhere, and plausibly why `veto_outcome_audit` reads 242/242 decided-against: the parity veto
cannot tell a harmless odd relabeling on an achiral complex from a real reflection, so it
declines them all. **The missing input is an achirality test at encode time — and `E` has the
3D coordinates to run one.** Designing that is not this chunk's job.

Open, not resolved here: 21/200 (10.5%) ruler-chiral passes whose string does *not* change under
reflection — an injectivity hole in `E`, or a ruler false-chiral near the `SPHERE_GAP` line.

## What this does to the 22.84-pt gap

| Block | pts | Roadmap says | Census C1 says |
|---|---|---|---|
| `slot_renumber` 252 | 5.04 | 201 enantiomers → v0.4.17 L1 CONSTRUCTION | **234 encoder (reflection non-canonical on achiral)**; 18 other |
| `structural` 484 | 9.68 | 301 DETACHED construction · 141 perception | 331 G-construction ✓ · **142 generator-correct** (E2/P) · 11 stereo |
| `rdkit_canonical` 113 | 2.26 | η-set drift | 62 G-construction · 47 generator-correct · 4 mirror |
| `facmer` 15 | 0.30 | arrangement | 13 G-construction (sphere), 2 generator-correct |
| `hard_fail` + `encode_fail` 278 | 5.56 | nothing produced | 252 no structure (unchanged; needs C3 parse-back) |

Generator-correct-but-failed, summed (`ok:same` + `ok:achiral` over the failing buckets):
236 + 142 + 47 + 2 + 12 (`hard_fail` with a good structure) = **439 molecules = 8.78 pts of the
22.84 are not the generator's fault.**

## Flagged passes (upper bound on false passes under the honest metric)

| Flag | n | Reading |
|---|---|---|
| MIRROR | 40 | all 40 are an inverted metal sphere with a byte-identical string → P1 (metal Δ/Λ) blind spot. Firm. Among passing chiral spheres: 1,643 same hand vs 112 inverted. |
| DIASTEREO | 104 | 18 sphere arrangement; 86 tetrahedral centre(s) differ. Unverified by eye. |
| G:sphere | 159 | 76 haptic, 83 σ; stored `coordination` calls 141 of them intact. ≥ 46 expected from the 1.2% squeeze floor alone. |
| G:ligand | 77 | likely deep clashes; unverified. |

Only the 40 should be quoted as false passes today. The rest is a work-list for C4.

## Side finding: UFF-tier bonds are 15% short

| | median M–L / Σr_cov | n |
|---|---|---|
| crystal (input) | 0.965 | 4,696 |
| generated, `UFF_1` | **0.818** | 4,643 |
| generated, `g-xTB_1` | 0.962 | 13 |

62.4% of generated structures have M–L bonds > 10% shorter than their crystal structure.
`KNOWN_LIMITATIONS.md` discusses *placement* targets (`bond_lengths.py`, `ENABLED_METALS`); this
is the geometry **after** UFF relaxation, which is what users receive. It does not hurt the
round trip (graph-ISO passes are *shorter*, 0.853 of crystal, than graph-ISO `structural`
failures, 0.990) — which is exactly why no sweep ever surfaced it. It is also what produces the
490 `ISO_CLASH` verdicts and the second-shell contacts that broke ruler version (a).

## Not done / next

- C2 (I2) should run the mirror mode over **all 5,000** — the 9.7% base rate above is a
  600-sample. It also owes renumbering and coordinate-noise stability, still unmeasured.
- The 142 generator-correct `structural` failures need the input-vs-generated *perceived*
  bond-order diff (plan rule 6) to split perceiver-fragility from geometry.
- Nothing here is frozen to `measurements/` yet; `results-census/` is gitignored. Freeze at C4.
