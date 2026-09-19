# The census — every molecule of the n=5,000 cohort attributed to one fault

**Date:** 2026-09-17 · **Branch:** `research/census` · **Plan:** [`CENSUS_PLAN.md`](CENSUS_PLAN.md)
**Chunks:** [C1 generator vs input](CENSUS_C1_generator_vs_input.md) ·
[C2 encoder self-consistency](CENSUS_C2_encoder_selfconsistency.md) ·
[C3 string sufficiency](CENSUS_C3_string_sufficiency.md) · C4 (this note, the join).
**Instrument:** `tools/census/attribution_table.py`. **Population:** the frozen v0.4.14 baseline
sweep, all 5,000 molecules, passes included. **Frozen:** `measurements/v0.4.17-census/`
(`attribution_table.tsv.gz` is the 5,000-row table; `census_*.jsonl.gz` the four instruments'
per-molecule verdicts). Every rate below is a count over that table; 1 molecule = 0.02 pts.

## What the census is

Four instruments, none of which routes through the encoder to judge the generator or through the
generator to judge the encoder:

| | judges | against | never imports |
|---|---|---|---|
| **I1** the neutral ruler (C1) | the generated XYZ | the input XYZ | `oinsmiles` |
| **I2** self-consistency (C2) | `E` | itself, under 10 transforms of the input | the generator |
| **I3** parse-back (C3) | `smiles_1` as a graph, no 3D | the input's ruler graph | the generator |
| **I4** the join (C4) | every molecule | the first rule that fires, in pipeline order | — |

The rule order is the pipeline's: no string → the string is wrong → nothing built → the wrong
graph built → the wrong hand built → the right structure re-read differently. A pass that trips
a rule is a **metric false pass**. Verification, printed by the tool before writing: 5,000 rows;
the bucket column reproduces `bucket_report_honest` exactly (3858/484/365/266/15/12); outcome
PASS 3,858 / FAIL 1,142; **`UNATTRIBUTED` = 0**; ruler `UNDECIDED` = 0.

## The partition of the 22.84-point gap

| fault | n | pts | mechanism (measured, not named after a bucket) | owner |
|---|---:|---:|---|---|
| **`G_CONSTRUCTION`** | **396** | **7.92** | generated heavy-atom graph ≠ input: 331 `DETACHED` (E2 wrote a slot-less fragment or the audit says detached), 60 donor set, 5 ligand bond | **generator** |
| **`E1_NONCANONICAL`** | **266** | **5.32** | structure and chirality reproduced; the input is achiral by the ruler and `E(mirror x) ≠ E(x)` — the parity veto fires on one hand (C2). 235 of these are `key_equal/slot_renumber` | **encoder** |
| `G_NOTHING` | 104 | 2.08 | no structure and a faithful string: 95 timeouts, 9 no-conformer | generator / budget |
| `E2_P_FRAGILE` | 81 | 1.62 | structure and chirality reproduced; the input's own string changes under renumbering or 0.02 Å noise (73 change the key) | perception |
| `E1_HCOUNT` | 80 | 1.60 | the string's H count ≠ input at the adapter level (bare non-donor atoms re-protonated on re-parse); 71 are `hard_fail` | encoder (serializer) |
| `E1_GRAPH` | 80 | 1.60 | the string's heavy graph ≠ input: 41 ligand bond, 39 donor set; 59 are `hard_fail` | encoder / perception |
| `E2_P_OTHER` | 67 | 1.34 | structure, chirality, H and self-consistency all fine, string still differs: bond orders / aromaticity / labels on a correct structure | perception |
| `P_DETACHED` | 18 | 0.36 | a ligand bound by distance, written as a slot-less fragment (oxo as free `[O-2]`, THF, DCM, η²-ethylene, the 3 `NON`) | perception |
| `P_E1_COVERAGE` | 12 | 0.24 | no string: encode timeouts | encoder |
| `G_DIASTEREOMER` | 11 | 0.22 | graph =, a different diastereomer (8 centres, 3 arrangement) | generator |
| `DATA_MULTI` | 10 | 0.20 | a second ion or molecule in the input file, unbound by distance | dataset |
| `G_HANDEDNESS` | 9 | 0.18 | graph =, the mirror image, and `E` does separate the hands | generator |
| `G_HCOUNT` | 7 | 0.14 | graph =, H decoration changed on the generated structure | generator |
| `E1_NONINJECTIVE` | 1 | 0.02 | the mirror image, and `E(mirror x) == E(x)` | encoder |
| **sum** | **1,142** | **22.84** | | |

**By owner:** generator **527 / 10.54** · encoder (E1) **427 / 8.54** · perception (E2) **148 /
2.96** · perception on the input (P) **30 / 0.60** · dataset **10 / 0.20**.

## What it refutes, in the roadmap's own terms

1. **"201 enantiomers → v0.4.17 L1 construction" is 9 molecules.** `G_HANDEDNESS` fails 9 times.
   The 252 `slot_renumber` are `E1_NONCANONICAL` 227 times: achiral molecules the encoder writes
   two ways. The lane that was sized at 4.02 pts of construction is a 5.32-pt **offline encoder
   lane** with one named mechanism, `oin/fold_parity.py::resolve` (C2).
2. **"`hard_fail` 266, 262 produce nothing ⇒ the throughput floor" is 104.** 158 of the 266 have a
   string that does not describe the input (`E1_HCOUNT` 71, `E1_GRAPH` 59, `P_DETACHED` 18,
   `DATA_MULTI` 10), and 84 of the 179 `hard_fail` timeouts are the generator working on the wrong
   molecule. The floor is 2.08 pts, not 5.24, and a serializer lane precedes it.
3. **`structural` is four owners, not one.** 313 generator (`DETACHED` dominates), 101 perception
   (`E2_P_FRAGILE` 65 + `E2_P_OTHER` 36), 30 encoder canonicality, 22 serializer, 18 other
   generator. The v0.4.16 estimate "INTACT+BOUNDARY = 82% perception, 141 molecules / 2.82 pts"
   is confirmed by an independent path: `E2_P_*` totals **148 / 2.96** over all buckets.
4. **`key_equal` is not "benign canonicalization — the win reclaimed".** 67 of its 365 members
   have a generated graph that differs from the input; the lossy key could not see it.
5. **`DETACHED` is 331, not 301** — the audit's DETACHED class among failures (367) overlaps
   `G_CONSTRUCTION` at 317; the census counts detachment from E2's own string (a fragment written
   with no slot) and the ruler's graph, which is stricter on both sides.

## The metric — 396 passes the round trip could not verify

| the pass hides | n | what happened |
|---|---:|---|
| `G_CONSTRUCTION` | 201 | the string was right, the generator built a different graph, and E2 read it back identically — `E` is not injective across graphs. ⚠ The ruler's own floor (C1 squeeze control 1.16%) says up to ~45 of these may be marginal-contact calls |
| `E1_NOT_ENCODED` | 106 | the generator changed a stereo element (99 a diastereomer, 7 the mirror) and the string never carried it. Ligand sp3 centres are not enforced **by design**; this is the notation's stated scope, counted so it is not mistaken for accuracy |
| `E1_GRAPH` + `E1_HCOUNT` | 56 | the string was wrong, the generator built it faithfully, E2 agreed |
| `E1_NONINJECTIVE` | 33 | the mirror image, byte-identical string (P1) |

**3,462 of 3,858 passes (69.24% of the cohort) are verified on the heavy-atom graph and on every
stereo element the ruler can see.** The 77.16% headline is a self-agreement rate. The gap between
the two, 7.92 pts, is not accuracy the project has; it is accuracy the metric asserts. The number
to carry forward is **69.24% verified, 77.16% self-consistent, 22.84 pts attributed**.

## What changes on the ladder

| release | was | is, by the census |
|---|---|---|
| **v0.4.17** | CONSTRUCTION: L1 201 enantiomers (4.02), L2 301 DETACHED (6.02) | **ENCODER CANONICALITY**: `E1_NONCANONICAL` 266 / 5.32 — a geometric achirality test in `fold_parity.resolve`, offline, mirror-audited against C2's 1,623 ruler-chiral inputs (36 protected pairs must stay separated). Then **L2 DETACHED 331 / 6.62** as construction. L1 handedness is 9 molecules and is dropped |
| **v0.4.18** | THE FLOOR: `hard_fail` 266 (5.32) + NO_STRUCTURE + facmer + encode_fail | **THE SERIALIZER, then the floor**: `E1_HCOUNT` 80 + `E1_GRAPH` 80 = 3.20 pts testable offline by parse-back (`tools/census/string_sufficiency.py parseback` is the gate), `P_DETACHED` 18 + `DATA_MULTI` 10 + `P_E1_COVERAGE` 12 = 0.80 as coverage; the compute floor is `G_NOTHING` 104 / 2.08 |
| **v0.4.19** | PERCEPTION 141 / 2.82 (estimate) | **PERCEPTION 148 / 2.96 (measured)**: `E2_P_FRAGILE` 81 (renumber/noise on the input — a P defect with no generator involved) + `E2_P_OTHER` 67 |
| the metric | `byte_exact` = the headline | **carry two numbers**: verified 69.24% and self-consistent 77.16%; every lever's sweep must report its effect on BOTH, or a lever that teaches E2 to agree with a wrong structure reads as a gain |

## Seed variance and the held-out 1,000 (C5, C6) — not run

The plan's two fresh-run probes need a quiet box and ~10 CPU-h. They are now stratifiable by
this table (`fault` column) and remain the way to put an error bar on 77.16% / 69.24%. Not
needed to choose v0.4.17: the top two classes are 12.5× and 8.4× larger than the seed error
their own sizes imply.

## Sources

`results-census/attribution_table.tsv` (5,000 rows, 37 columns — one per instrument flag, so the
overlaps the first-match rule hides can be read back), `attribution_summary.json`; frozen under
`measurements/v0.4.17-census/`. The instruments: `tools/census/{neutral_graph, g_vs_input,
mirror_probe, e_selfconsistency, veto_probe, string_sufficiency, attribution_table}.py`.
