# v0.4.18 — the census, re-run on the v0.4.17 sweep

> **17.28 pts self-consistent / 25.40 pts VERIFIED, every molecule attributed, `UNATTRIBUTED 0`.**
> **E1 fell 8.54 → 3.38. G is now 56% of what fails. DETACHED is 343 molecules / 6.86 pts and 86% of
> it is η-bound.** L1d — the residual non-canonical strings — is worth **0.24 pts** and is demoted.

Branch `research/v0417-reattribution`, 2026-09-19, off `main` @ `v0.4.17`. Owner request: re-attribute
before choosing the next lane. The census table was stale for every molecule `OIN_EXACT_DONOR_FOLD`
touched, and both gap numbers leaned on it.

## The partition

FAIL = 864 molecules = **17.28 pts** (was 1,142 / 22.84):

| owner | fault | old | new | pts | what it is |
|---|---|---:|---:|---:|---|
| **G** | `G_CONSTRUCTION` | 396 | **386** | **7.72** | generated graph ≠ input. **`DETACHED` 334**, `SPHERE_DIFF` 49, `LIGAND_DIFF` 3 |
| | `G_NOTHING` | 104 | **76** | 1.52 | no structure: `TIMEOUT` 68, `NOCONF` 8 |
| | `G_HCOUNT` · `G_DIASTEREOMER` · `G_HANDEDNESS` | 27 | 19 | 0.38 | 7 · 7 · 5 |
| | **G total** | **527** | **481** | **9.62** | |
| **E2** | `E2_P_FRAGILE` | 81 | **101** | 2.02 | string changes under renumbering / 0.02 Å noise (93 at key level) |
| | `E2_P_OTHER` | 67 | **74** | 1.48 | bond orders / labels on a correct structure |
| | **E2 total** | **148** | **175** | **3.50** | |
| **E1** | `E1_HCOUNT` | 80 | 79 | 1.58 | the string's H count ≠ the input's |
| | `E1_GRAPH` | 80 | 77 | 1.54 | the string's graph ≠ the input's |
| | `E1_NONCANONICAL` | **266** | **12** | 0.24 | one molecule, two strings |
| | `E1_NONINJECTIVE` | 1 | 1 | 0.02 | |
| | **E1 total** | **427** | **169** | **3.38** | |
| **P** | `P_DETACHED` 18 · `P_E1_COVERAGE` 11 | 30 | 29 | 0.58 | |
| **data** | `DATA_MULTI` | 10 | 10 | 0.20 | two molecules in one file |

VERIFIED gap = FAIL + the **406 passes that carry a fault** (was 396) = **1,270 molecules = 25.40 pts**:

| false pass | old | new | |
|---|---:|---:|---|
| `G_CONSTRUCTION` | 201 | 200 | the round trip agrees about a structure whose graph differs: `SPHERE_DIFF` 128, `LIGAND_DIFF` 63, `DETACHED` 9 |
| `E1_NOT_ENCODED` | 106 | 109 | G changed a stereo element the string never carried |
| `E1_GRAPH` · `E1_HCOUNT` | 56 | 60 | the string was wrong and G built what it said |
| `E1_NONINJECTIVE` | 33 | 37 | mirror pair, one string |

By owner, verified: **G 681 / 13.62 · E1 375 / 7.50 · E2 175 / 3.50 · P 29 / 0.58 · data 10 / 0.20.**

## What it says about the ladder

1. **L1d is demoted.** `E(mirror x) ≠ E(x)` on a ruler-achiral molecule still holds for **120**
   (was 670). 75 of them pass anyway, 33 fail for a reason that fires earlier, and **12 fail because
   of it — 0.24 pts.** It remains a *canonicality* defect (one molecule, two hashes) and should be
   fixed for that reason, with an automorphism-based achirality bit; it is not an accuracy lane.
2. **L2 DETACHED is confirmed and sharpened: 343 molecules / 6.86 pts** (334 fail, 9 pass).
   **295 of 343 are η-bound (86%)**; geometry `TET` 163, `TPL` 55, `OCT` 34, `SPY` 28, `TPY` 28 —
   piano stools. It is 39% of everything that fails and the largest single block by a factor of 3.
3. **The serializer (v0.4.18 on the old ladder): `E1_HCOUNT` + `E1_GRAPH` = 156 fail / 3.12 pts,
   216 / 4.32 verified.** 128 of the 216 sit in `hard_fail` — a generator asked to build the wrong
   molecule — so part of "the compute floor" is this, not throughput.
4. **E2 (perception) grew, 2.96 → 3.50, and nothing got worse.** `E1_NONCANONICAL` fires before the
   E2 rules; 25 molecules the lever did not rescue were E2 all along, hidden behind it.
5. **The owner's parked lane, sized:** the generator produced nothing, or the encoder timed out, for
   **87 molecules = 1.74 pts** (`G_NOTHING` 76 + `P_E1_COVERAGE` 11). That is the ceiling on what
   an always-works fallback can add to self-consistent accuracy — and a fallback structure still
   has to get past the neutral ruler to count as verified.
6. **Stereo the string does not carry: 146 passes / 2.92 verified pts** (`E1_NOT_ENCODED` 109 +
   `E1_NONINJECTIVE` 37). Invisible to the self-consistent number. Y1 P3 lives here.

**Recommended next lane: L2 DETACHED, scoped to η.** Then the serializer.

## Where the faults went

| | molecules | keep their fault | change |
|---|---:|---:|---:|
| movers (the lever's complete 504) | 504 | 188 | **316** |
| non-movers | 4,496 | **4,482** | **14 (0.3%)** |

Movers: **223 `E1_NONCANONICAL → NONE`** (not merely passing — *verified*), 19 `G_CONSTRUCTION → NONE`,
14 `G_NOTHING → NONE`, 19 `E1_NONCANONICAL → E2_*`, 4 `→ G_CONSTRUCTION`. `PASS → FAIL` 3,
`FAIL → PASS` 279.

**The 14 non-movers are the noise floor of attribution itself: 0.28 pts.** 6 of them are
`E1_NONCANONICAL → E2_*`, which is the lever after all — it moved their *mirror's* string, not
theirs, so they are non-movers by the A/B's definition and movers by rule 6's. The other 8 are
timing: `G_NOTHING → {E2_P_FRAGILE 4, G_CONSTRUCTION 1, NONE 1}`, `P_E1_COVERAGE → G_NOTHING` 1,
`E1_NOT_ENCODED → E1_NONINJECTIVE` 1. **Do not choose a lane smaller than ~0.3 pts from this table.**

## Instruments, and what a broken one would print

`tools/census/attribution_table.py` was hard-wired to the v0.4.14 sweep. It now takes the sweep and
one path per instrument, and:

- **a missing instrument file is FATAL.** It used to read as zero rows, and a join over zero rows
  is a plausible table with one class quietly absent.
- **the bucket column must reproduce that sweep's frozen `bucket_report_honest`**, and
- **`NONE` must equal the VERIFIED count `tools/v0417/sweep_two_numbers.py` printed for that sweep**
  (a six-line re-implementation of rules 1–5). Two tools, one answer: 3,462 and **3,730**. Both agree.

**CONTROL, run first:** the modified tool on the *record* sweep reproduces the census table
**byte-for-byte** — against `results-census/attribution_table.tsv`, against the **public frozen**
`measurements/v0.4.17-census/attribution_table.tsv.gz`, and the summary JSON.

Inputs, by whether they depend on the sweep:

| instrument | source | |
|---|---|---|
| honest buckets, C1 `g_verdict`, C3 `parseback` | `results-v0.4.17-sweep/` | from `post_sweep.sh` |
| C2 `e_selfconsistency` | `results-v0.4.17-exactfold/e_selfconsistency_exact.jsonl` | the live full-cohort run under the exact fold = the shipped encoder |
| C3 `parseback_gen` (4,780 structures), `pflags --charge-probe` (5,000; probe control 1,035/1,035), attach-class audit (5,000 reports, 0 missing) | **derived here** | |
| C1 mirrored-input ruler verdict, natural-twin clusters | `results-census/` | properties of the INPUT — reusable |

`tools/census/reattribution_diff.py` is new: old table vs new, split by mover. The split *is* the
control — a wrong join key shows up as thousands of non-movers changing fault. 14 did.

## Not done

- C5 seed variance and C6 held-out cohort remain optional, as they were.
- `pflags` hit its 120 s per-molecule limit on 19 inputs. It feeds descriptive columns only; no rule
  reads it.
- The 200 `G_CONSTRUCTION` false passes were not split further. `SPHERE_DIFF` 128 is the ruler
  disagreeing about the coordination sphere on a structure the round trip accepts, and whether
  that is G or the ruler's tolerance is a question for the L2 lane.
