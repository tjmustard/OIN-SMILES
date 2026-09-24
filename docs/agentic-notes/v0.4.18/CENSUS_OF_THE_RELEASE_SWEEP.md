# v0.4.18 — the census, re-run on the release sweep

> **13.06 pts self-consistent / 21.60 pts VERIFIED, every molecule attributed, `UNATTRIBUTED 0`.**
> **The noise floor of attribution is zero: all 4,157 non-movers keep their fault.** G fell 481 → 278;
> `DETACHED` 334 → 142. The three largest blocks are now within 1 pt of each other: G construction
> 3.94 · E2 perception 3.52 · E1 serializer 2.96.

Branch `research/v0418-eta-selection`, 2026-09-23, under the owner's delegation. Run by
`tools/v0418/run_census_release.sh` on `results-v0.4.18-release-sweep` (86.94% / 78.40%). The
census must describe what *ships* before the next lane is chosen; both v0.4.18 promotions moved
what the generator builds, and the v0.4.17 table was stale for 843 molecules.

## Controls, first

- The parameterised tool on the v0.4.14 record sweep reproduces the frozen census table
  **byte-for-byte** (`cmp` against `measurements/v0.4.17-census/attribution_table.tsv.gz`).
- `NONE` = 3,920 = the VERIFIED count `sweep_two_numbers.py` printed for this sweep. Two tools, one answer.
- **Movers are defined by coordinates**: a molecule whose generated structure differs from the v0.4.17
  sweep's (843). A non-mover gets the same structure, the same two strings, and — this time — the same
  fault: **4,157 / 4,157**. The v0.4.17 re-attribution had 14 non-movers change (0.28 pts of noise);
  here there are none, because the release sweep is byte-identical to the v0.4.17 sweep outside the
  molecules the levers touch.

Inputs: C1 ruler + C3 parse-back of `smiles_1` from `post_sweep.sh`; C3 parse-back of the generated
side, perception flags with the charge probe, and the attach-class audit derived here; C2 encoder
self-consistency **reused** from the v0.4.17 exact-fold run (the encoder is unchanged — `smiles_1`
moved on 0 of 5,000 in every sweep since); the mirror ruler and twin clusters from the census.

## The partition

FAIL = 653 molecules = **13.06 pts** (v0.4.17: 864 / 17.28):

| owner | fault | v0.4.17 | **v0.4.18** | pts | what it is |
|---|---|---:|---:|---:|---|
| **G** | `G_CONSTRUCTION` | 386 | **197** | **3.94** | **`DETACHED` 142** (151 with its 9 false passes: η 103 / non-η 48), `SPHERE_DIFF` 50, `LIGAND_DIFF` 5 |
| | `G_NOTHING` | 76 | 69 | 1.38 | `TIMEOUT` 61, `NOCONF` 8 |
| | `G_HCOUNT` · `G_DIASTEREOMER` · `G_HANDEDNESS` | 19 | 12 | 0.24 | 3 · 5 · 4 |
| | **G total** | **481** | **278** | **5.56** | |
| **E2** | `E2_P_FRAGILE` | 101 | 105 | 2.10 | string changes under renumbering / 0.02 Å noise (`key` 96) |
| | `E2_P_OTHER` | 74 | 71 | 1.42 | `NO_STEREO` 42, `SAME` 29 |
| | **E2 total** | **175** | **176** | **3.52** | |
| **E1** | `E1_HCOUNT` | 79 | 77 | 1.54 | `h_diff=1` 40 · 2 10 · 3 13 · 4 6; 72 of 84 sit in `hard_fail` |
| | `E1_GRAPH` | 77 | 71 | 1.42 | `LIGAND_DIFF` 39, `SPHERE_DIFF` 32 |
| | `E1_NONCANONICAL` | 12 | 11 | 0.22 | |
| | **E1 total** | **169** | **160** | **3.20** | |
| **P** | `P_DETACHED` 18 · `P_E1_COVERAGE` 11 | 29 | 29 | 0.58 | |
| **data** | `DATA_MULTI` | 10 | 10 | 0.20 | |

VERIFIED gap = FAIL + **427 passes that carry a fault** (was 406) = 1,080 molecules = **21.60 pts**:
`G_CONSTRUCTION` 205 (`SPHERE_DIFF` 133, `LIGAND_DIFF` 63, `DETACHED` 9), `E1_NOT_ENCODED` 118,
`E1_GRAPH` 61, `E1_NONINJECTIVE` 36, `E1_HCOUNT` 7. By owner, verified: **G 483 / 9.66 · E1 382 / 7.64
· E2 176 / 3.52 · P 29 / 0.58 · data 10 / 0.20.**

## Where the faults went (movers 843, non-movers 4,157)

Movers: 559 keep their fault, 284 change — **177 `G_CONSTRUCTION → NONE`** (verified, not merely
passing), 19 `E2_P_OTHER → NONE`, 7 `G_NOTHING → NONE`, 3 each `G_HCOUNT` / `E1_NONCANONICAL` /
`E1_NONINJECTIVE → NONE`; the other way 18 `NONE → G_CONSTRUCTION`, 4 `→ E2_P_OTHER`, 3 `→ E1_NOT_ENCODED`.
13 `G_CONSTRUCTION → E2_P_OTHER` and 8 `→ E1_NOT_ENCODED`: a ring that now attaches exposes the fault
behind it. `PASS → FAIL` 16, `FAIL → PASS` 227. Non-movers: 0 change.

## What it says about the ladder

1. **The residual `DETACHED` (142 fail, 151 with false passes) is two populations.** η-bound 103 —
   Ru 37, Ti 12, Zr 11, Rh 9, Fe 8 — mostly the ~30 molecules of oracle headroom the per-atom arrival
   count could not decide, plus rings the retarget did not reach; and **non-η 48** (Y 7, Pt 6, Ir 4,
   Pd 4, Zn 4), which no v0.4.18 lever touched: a σ-donor that never arrives. `TET` 72 of 151.
2. **`SPHERE_DIFF` is the largest *verified* block G owns: 183 (133 of them false passes).** The ruler
   disagrees about the coordination sphere on a structure the round trip accepts. Whether that is G
   or the ruler's tolerance is still the open question the v0.4.17 census left.
3. **E2 perception (176 / 3.52) is now within 0.4 pts of G construction.** `E2_P_FRAGILE` 105 is 92
   non-η, Ni 24 / Zn 22: strings that change under renumbering or 0.02 Å noise — an encoder
   *stability* defect, one molecule, two hashes.
4. **The serializer (`E1_HCOUNT` + `E1_GRAPH` = 148 fail / 2.96 pts; 216 verified / 4.32)** is mapped
   in `spec/handoffs/v0.4.19/SERIALIZER_LANE_MAP.md`: three different defects, only one of them in the
   string-writing code (an η-rank swap in the aligner; 5 of 12 sampled rows flip to ISO when its
   fail-safe is forced), the others charge neutralisation writing a 0-H anion bare and perception's
   valence cap deleting ligand bonds. 129 of the 148 sit in `hard_fail`.
5. **`E1_NOT_ENCODED` grew 109 → 118**: stereo the string does not carry, invisible to the
   self-consistent number. Y1 P3 lives here (2.36 verified pts with `E1_NONINJECTIVE`).
6. The owner's parked lane: `G_NOTHING` 69 + `P_E1_COVERAGE` 11 = **80 molecules / 1.60 pts** ceiling.

**Recommended next lane (v0.4.19), for the owner:** E2 perception fragility (105 / 2.10, non-η,
Ni/Zn) or the serializer's η-rank swap — both are *encoder* lanes and neither needs a generator sweep
to measure (the offline instruments are exact for string-side changes). The residual DETACHED needs
a different selector than the arrival count and is smaller per molecule of effort now.

## Instruments

`tools/v0418/run_census_release.sh` (control → attach → parse-back gen → pflags → table → movers →
diff). Frozen: `measurements/v0.4.18-release-sweep/` (`v0418_rcensus_*`).
