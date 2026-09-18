# Census C3 — is the string sufficient? Parse-back, natural twins, perception flags

**Date:** 2026-09-17 · **Branch:** `research/census` · **Plan:** [`CENSUS_PLAN.md`](CENSUS_PLAN.md) ·
**Before this:** [`CENSUS_C1_generator_vs_input.md`](CENSUS_C1_generator_vs_input.md),
[`CENSUS_C2_encoder_selfconsistency.md`](CENSUS_C2_encoder_selfconsistency.md)
**Instrument:** `tools/census/string_sufficiency.py` (four subcommands: `parseback`, `collide`,
`pflags`, `summary`). No generator anywhere. Parse-back and the collision scan each take ~15 s on
4 workers; the perception pass with the charge probe takes ~1 h on 6 workers (5,000 perceptions +
2 × 1,038 full encodes).
**Population:** all 5,000 cohort inputs (4,988 with a string) and, for the generated-side control,
the 4,748 that have a generated structure in the v0.4.14 sweep.
**Outputs (gitignored, freeze at C4):** `results-census/{parseback, parseback_gen,
parseback_control-delete, parseback_control-reattach, pflags}.jsonl`, `collisions.json`,
`string_sufficiency_summary.json`, and the `.log` beside each.

## Headline

1. **The string describes the input for 97.0% of encodable molecules (4,838/4,988) and for
   only 40.6% of `hard_fail` (108/266).** Reading `smiles_1` back into a graph with no 3D and
   comparing it to the input's neutral graph (the I1 ruler's, `oinsmiles` never imported on that
   side): 150 strings carry a different heavy-atom graph, 28 carry a fragment with no slot that
   the adapter refuses, and 84 carry an H count the generator cannot recover. Those three faults
   are almost entirely a `hard_fail` phenomenon:

   | `hard_fail` (n=266) | n | % | what the sweep recorded |
   |---|---|---|---|
   | string faithful (graph ISO, H same, all fragments slotted) | 108 | 40.6 | 95 timeouts, 9 generation failures, 4 blank |
   | **H count differs at the adapter level** | **71** | 26.7 | 51 timeouts, **16 "Atom count mismatch"**, 4 gen. failures |
   | ligand graph differs (a covalent bond made or broken) | 37 | 13.9 | 19 gen. failures, 17 timeouts |
   | **a fragment with no slot** (adapter raises before embedding) | **28** | 10.5 | 28 gen. failures |
   | metal sphere differs (donor set) | 22 | 8.3 | 16 timeouts, 5 atom-count mismatches |

   **59.4% of `hard_fail` (158/266) — 3.16 pts of the 22.84-pt gap — has a string that does not
   describe the input**, and 84 of the 179 `hard_fail` timeouts (47%) are the generator working
   on the wrong molecule, not a throughput problem. The roadmap's "266 produce nothing ⇒ v0.4.18
   floor lane" is therefore mostly a **serializer (E1) lane**, with a perception sub-split below.
2. **The 84 H-count faults have one mechanism, already in project memory as a hypothesis:** a
   non-donor atom that perception left with an unusual valence (a `[CH]` radical ring, a bare `C`
   with three heavy neighbours in a non-aromatic ring) serializes bare and is silently
   re-protonated when the adapter re-reads it — the adapter's donor-H reconciliation repairs the
   *donor* convention (parser-level H differs for 3,278 molecules, adapter-level for 84) but has
   no rule for non-donors. 71/84 are `hard_fail` (85% conditional failure), +1 to +14 H each.
3. **The 28 slot-less fragments split 10 / 18.** Ten are genuinely separate molecules in the
   input file — `[BH4]⁻` ×4, `BF4⁻`, triflate, acetate, a free quinoxaline ligand ×2, a lone
   `[H]` — the ruler agrees they are unbound (graph ISO, same component count): a **dataset /
   coverage** fault, the OIN format has no way to write an uncoordinated ion. Eighteen are bound
   by distance and **detached by perception**: an oxo written as free `[O⁻²]` ×5, a coordinated
   THF, a coordinated DCM, an η²-ethylene, `B(C₆F₅)₃` adducts ×2, a carborane, and the three
   `NON`-geometry strings whose metal has no donor at all (CN 10 by distance, 0 in the string).
4. **The parse-back settles C1's 236 `byte_exact` passes whose generated graph was not the
   input's: 201 have a correct string and 35 do not.** For the 201, E1 was right, G built a
   different graph and E2 read it back identically — **E is not injective across graphs and the
   metric cannot see it (5.2% of passes)**. For the 35, E1 wrote the wrong graph, G built it
   faithfully and E2 agreed: a consistent wrong answer. A further 17 `byte_exact` strings read
   non-ISO while C1 says G reproduced the input — that is the ruler's own marginal-contact floor
   (0.3%), and it bounds how much of the 150 is instrument rather than serializer.
5. **In the wild `E` never collided and was non-canonical for 3 of 47 natural-twin clusters.**
   The 1,033 cat/photo twins are byte-identical files (checked; E is deterministic by C2), so the
   scan clusters the 5,000 inputs by heavy-graph isomorphism instead: 47 clusters, 110 molecules,
   102 pairs judged by the ruler input-vs-input. **0 pairs share a string across a ruler
   MIRROR/DIFFERENT verdict** (no collision, including the 18 `_comp_n` sibling refcodes, 3 of
   which are enantiomer pairs the string separates). **50 pairs are one molecule with two strings**,
   45 of them one cluster: **eleven crystal structures of Zn tetraphenylporphyrin encode to six
   different strings** — different Kekulé/aromatic assignments of the macrocycle, phenyls written
   as `[CH]` radical rings — and all eleven fail the round trip. The other two clusters are a
   dppz-type diimine perceived aromatic in one structure (fails) and quinoid in another (passes),
   and the C1/C2 slot-label swap on a symmetric bis-phosphine (`P{2}…P{3}` vs `P{3}…P{2}`, one
   passes, one `key_equal`). Bond-order perception and slot labelling, the two mechanisms C2
   named, are the two seen across independent structures.
6. **The encoder ignores the total charge.** `XYZToSMILES.convert` sets `charge = 0` and never
   reads the comment line's `Charge:` field; perception then assigns the metal `0 − Σ(ligand
   charges)`. 1,038 cohort inputs (20.8%) state a non-zero charge and 2,078 state none. The metal's
   formal charge feeds `metal_d_electron_count` → the electronic geometry-template override, so
   the wrong charge *can* change a geometry code. **It almost never does: honouring the stated
   charge changes 10 of 1,036 strings** (all geometry-code flips, none a ligand charge), while
   the perceived oxidation state is implausible for 70% of charged inputs. A lossiness worth one
   line in `translator.py`, not a lane — **§ Perception flags** below.

## What the instrument does, and what a broken one would print

**Parse-back.** Two graphs are built from `smiles_1` and compared to the input's ruler graph by
heavy-atom isomorphism, with the ruler's own ladder on the XYZ side (metal reading std/min/max,
then the `cov=0.88` clash rescue; first match wins and is labelled `ISO` / `ISO_MARGINAL` /
`ISO_CLASH`). *Parser level* is the inline handler alone — fragments and `{n}` slot markers,
`MolFromSmiles(sanitize=False)`, nothing else. *Adapter level* is the ligand SMILES the generator
hands to MetalloGen after `_prepare_ligand_fragments` (kekulisation, the Cp charge rescue, the
donor-H reconciliation), donors read from the atom-map numbers. H counts are compared per
symmetry class as the ruler does and reported beside the graph verdict, never folded into it.
Hydride fragments `[H]{n}` fold into the metal's H count exactly as the ruler folds a
metal-bound H.

Controls, full cohort, read before any number above:

| control | must read | usable | ISO | verdict |
|---|---|---|---|---|
| `delete` — drop the last slotted heavy ligand | non-ISO | 4,985 | **0** | PASS |
| `reattach` — move ONE metal edge to a neighbour of its donor that is not symmetry-equivalent | non-ISO | 4,972 | 0 at std, **1** at the `max` metal reading (`XEWSID`, a Y–C contact inside the 1.10× cutoff) | PASS at the reading the real numbers use |

Both controls first *failed*, and each failure was the instrument telling the truth: `delete`
read ISO on 79 molecules whose last slotted fragment was a hydride (a heavy-atom graph cannot
see an `[H]` leave; all 79 read H-DIFF), and `reattach` read ISO on 11 where the neighbour is an
automorphic image of the donor in a bond-order-free graph — an unsubstituted phenyl, N₂, a
tetrazolate, the N3/N4 of a 1,2,4-triazolyl. Both controls now exclude what the graph is blind
to by construction, and the blind spots are stated here rather than absorbed.

A broken parse-back that ignored slot markers would still read ISO on every ligand graph and
would fail `reattach`; one that ignored fragments would fail `delete`. A parse-back whose atom
indices were off by one would read `SPHERE_DIFF` on most of the cohort, not on 79.

**Generated-side control.** The same parse-back run on `smiles_2` against the *generated* XYZ
(4,748 molecules): ISO-family 4,249 (89.6%), `SPHERE_DIFF` 417, `LIGAND_DIFF` 77, and **318
strings with a slot-less fragment (313 `structural`)** — E2 wrote the fragment with no slot
because the generator left it unbound. 289 of those 318 are the attach-class audit's DETACHED
molecules among genuine failures (315) — the cross-check the plan asked for, from an
instrument that never read the audit. Among `byte_exact` passes 281 (7.3%) strings do not match
the generated structure's distance graph although E2 re-read them identically: the generator
side of the metric's blindness, the same population as headline 4 seen from the other end.

**Collision scan.** Inputs are grouped by heavy formula with H, pairwise `isomorphism` inside a
group (123 tests), union-find into clusters; inside each multi-member cluster every member is
compared to the first with the full ruler (input vs input, so `MIRROR`/`DIFFERENT`/`SAME` mean
what they meant in C1) and every pair with differing strings is compared directly. Same string
across two clusters would be a graph collision: 0. The ruler's blind spots (E/Z, atropisomers,
planar chirality) can only turn a real difference into a false `SAME`, i.e. inflate
NON_CANONICAL — the three clusters above were read by eye and the string diffs are perception
and slot labels, not a stereo element the ruler missed.

## The non-ISO 150, by mechanism

`SPHERE_DIFF` (79): in 55 the input's distance graph has metal–C contacts the string lacks (one
C: 32; two: 10; four: 2) — agostic/close contacts and partial η-rings that perception did not
make bonds, 24 of them `byte_exact` (E2 read the generated structure the same way, so they cost
nothing). In 12 the string has a donor (O, N, C) the ruler does not see even at the 1.10×
cutoff — perception's metal cutoff is longer than the ruler's. `LIGAND_DIFF` (71): 43
`hard_fail`, 15 `byte_exact`; the `byte_exact` ones split a ligand into two components (7) or
merge (3) and pass because E2 makes the same call on the generated geometry. Marginal heavy-atom
contacts (within 6% of the covalent cutoff) are present in 31% of the non-ISO inputs vs 18% of
the ISO ones, and the ruler's second-shell filter drops 5× more contacts in them: the 150 are
enriched in geometries where distance and perception legitimately disagree, and the 17 of
headline 4 are the size of that floor.

## The "ligand as written fails RDKit sanitize" flag

21.7% of `byte_exact` strings (838/3,858) and 51.5% of `hard_fail` (137/266) contain a ligand
fragment RDKit refuses to sanitize as written. **Not usable as a fault flag on its own:** the
`KekulizeException` half is dominated by the neutral-Cp convention (`[cH]1[cH][cH][cH][cH]1`,
which the adapter charges and kekulizes deliberately) and the `AtomValenceException` half by
uncharged borates and P(V) ylides that the generator handles. Kept in the record for C4 as a
covariate, with this caveat.

## Perception flags — `pflags --charge-probe`

Every input perceived the way the encoder perceives it (`get_tmc_mol(path, 0)`): 4,985 OK, 15
perception timeouts at 120 s. Then, for each of the 1,036 inputs stating a non-zero charge, two
full encodes through a patched `perception_tmc.get_tmc_mol`: one at charge 0 (must reproduce
`smiles_1` — **1,036/1,036**, the patch path is the production path) and one at the stated
charge (the perceived total must equal it — **1,031/1,036**, the patch fires).

**The hardcoded charge is a real lossiness and a negligible round-trip lever.**

| | n | perceived metal charge in `ALLOWED_OXIDATION_STATES` |
|---|---|---|
| stated `Charge: 0` | 1,872 | **91.5%** |
| stated non-zero | 1,036 | **30.3%** as perceived · **71.4%** if the charge were honoured |
| no charge stated | 2,092 | 82.4% |

For the charged inputs the metal is off by exactly the stated charge — Ru(I) for 114 Ru(II)
cations, Ir(II) for 61 Ir(III), Cu(0)/Au(0)/Rh(0) for 58/43/55 monocations, Pd(I) for 55 — and
`ALLOWED_OXIDATION_STATES` would have caught 722 of them had anything called it. The remaining
296 charged inputs whose honoured charge is *still* implausible (Pd(0) 43, Pd(−II) 29, from
stated −2/−4) are cases where perception's per-ligand charge assignment, not the total, is the
suspect.

**But the total charge reaches only the metal's formal charge, and that reaches only the
d-electron geometry prior.** Honouring it changes the string for **10 of 1,036** (5 `byte_exact`,
3 `key_equal`, 2 `structural`, 0 `hard_fail`), every one a geometry-code flip (TPY→TET 5,
TBP→SPY 4, SPL→TET 1: Cu(I) and Rh(I) d¹⁰/d⁸ priors on 4- and 5-coordinate spheres) and never
a ligand charge — the per-ligand charges come from orbital filling and do not see the total.
At most ±5 molecules (±0.1 pt) either way. **Not a v0.4.17 lane.** The correct fix is one line
in `translator.py` for lossless *chemistry*; it needs a mirror audit and a sweep like any lever
but it will not move the headline.

**Perception plausibility as a covariate for C4** (perceived under charge 0, all 4,985):

| flag | prevalence | `byte_exact` | `structural` | `key_equal` | `hard_fail` |
|---|---|---|---|---|---|
| metal charge not in allowed set | 25.1% | 22.4% | 32.1% | 34.0% | **39.5%** |
| a ligand atom with a radical electron | 16.3% | 17.2% | 12.4% | 4.4% | **25.9%** |
| a charged carbon (X-type C donor) | 43.7% | 40.2% | **61.8%** | 47.1% | 54.5% |
| any atom with \|q\| ≥ 2 | 5.4% | — | — | — | — |

The first row conflates the hardcode with genuine mis-perception; restrict to the stated-neutral
1,872 (8.5% implausible) before using it as a fault flag. None of these is a fault on its own —
they are the ligand classes perception handles worst, to be crossed with C2's fragility (which
is *not* charge-correlated: 21.5% of charged inputs are renumber/noise-fragile vs 18.1% of the
rest). Pass rate by stated-charge class: `Charge: 0` 82.4%, non-zero 78.0%, none stated 72.1%
— the unstated group is the dataset's hard cases, not the encoder's.

## What this changes in the plan

- **Fault-table row 2 (E1 serializer/parser) is now measurable and is the largest single
  `hard_fail` cause**: graph ≠ input 59, slot-less fragment 28, H count ≠ input 71 — 158 of
  266. Row 3 ("no generated structure → G-nothing, timeout vs no-conformer") must be applied
  *after* row 2: 47% of the `hard_fail` timeouts are G given the wrong molecule.
- **Row 2 needs three sub-rows** for C4: `E1_GRAPH` (ruler-vs-string heavy graph, minus the ~17
  floor), `E1_HCOUNT` (adapter-level H decoration DIFF), `P_DETACHED` / `DATA_MULTI` (slot-less
  fragment, bound vs unbound by distance).
- **Row 7 (metric false pass) has two sides**, both counted now: G-side 201 (string right,
  structure wrong, E2 blind), E1-side 35 (string wrong, structure wrong the same way).
- The natural-twin scan supports C2's mechanism list and adds nothing new: perception bond
  orders and slot labels are the two ways one molecule gets two strings. No injectivity
  failure was found in the 102 pairs the cohort offers — the P1 blind spot (metal Δ/Λ
  enantiomers under a byte-identical string) is real (C1: 40 passes) but the twins do not
  sample it.
- The charge finding (headline 6) is **not a lane**: the probe bounds it at ±5 molecules. Fix
  it for lossless chemistry (one line in `translator.py`, plus wiring `ALLOWED_OXIDATION_STATES`
  into perception as a warning), mirror-audit and sweep as for any lever, and expect the
  headline not to move. The per-molecule flags (`os_in_allowed`, radicals, charged carbon)
  join C4 as covariates for the row-6 P sub-split, not as faults.

## Traps this chunk added or confirmed

- **A control can fail because it is right.** Both parse-back controls failed on the full cohort
  after passing on a 60-molecule sample: the graph is blind to hydrides and to automorphic
  re-attachment *by design*. Run controls on the full population, and when one fails, name the
  blind spot before touching the instrument.
- **`XYZToSMILES.convert` hardcodes `charge = 0`.** Every accuracy figure so far was measured
  with the total charge ignored; the string carries no metal charge either, so the round trip
  is self-consistent and blind to it.
- **`ALLOWED_OXIDATION_STATES` is defined in `perception_tmc.py` and used nowhere in `src/`.**
- The attach-class audit's DETACHED class includes 52 `byte_exact` passes (its own stated
  control) — use E2's slot-less fragment as the strict definition.
- `#DONE` trailers: every JSONL in `results-census/` ends with one; a reader that does not skip
  `#` lines dies on the last line.
- Never split a string on `.` and call the last slotted piece "a ligand" without checking it is
  a heavy atom.

## Next (C4)

I4 join: `tools/census/attribution_table.py`, one TSV of 5,000 rows joining the bucket,
`g_verdict` (C1), `e_selfconsistency` (C2), `parseback` + `pflags` (C3) and the attach class;
fault by first match with row 2 split as above; `CENSUS.md`; `/freeze-measurements` with the
harvester allowlist extended for `results-census/`; re-point `ROADMAP_100_100.md`.
