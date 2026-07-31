# v0.4.15 Lane 1 — `OIN_ATTACH_RETURN` · **REFUTED as an accuracy lever**

**Verdict: default OFF. It changes the returned structure 51 times across 289 molecules and
improves it zero times, at essentially no cost.**

## What it does

v0.4.7 built a coordinate-only attachment predicate (`conformer_ligands_attached`) and wired it
into **acceptance** only. `OIN_ATTACH_CHECK`'s own lever entry names the gap that left as its
third residual class:

> GAVSED (acceptance rejected everything, and `_select_by_geometry`'s fallback ranking never
> consults this check — **the check guards ACCEPTANCE, not RETURN**).

Lane 1 is that check on the return path. `_attach_rank` is a sort key demoting a detached
conformer below an attached one at **all three** exits of `_select_by_geometry_impl` — the
geometry-classified sort, the eta-winding scan, and the lowest-energy fallback — not just the one
the charter named. When no conformer holds its sites, the pre-lever answer is still returned, so
the lever cannot lower accuracy; `OIN_ATTACH_RETURN_STRICT` is a separate lever for converting
that case into an honest `hard_fail`.

## Measured

Pre-flight (`tools/attach_return_preflight.py`, `measurements/v0.4.15/`):

| population | n | `SITE_LOST` |
|---|---:|---:|
| `structural`/`DETACHED` | 301 | **289 (96.0%)** |
| `byte_exact`/`DETACHED` | 52 | **1 (1.9%)** |
| `byte_exact`/`INTACT` | 250 | **0 (0.0%)** |

A/B, real generation, honest scoring (`measurements/v0.4.15/lane1_*.json`):

| arm | n | gains | losses | moved | runtime | `>30 s` |
|---|---:|---:|---:|---:|---:|---:|
| target `SITE_LOST` | 289 | **0** | 0 | **51** | 1.01× | 147→148 |
| control `byte_exact`/DETACHED | 52 | 0 | **0** | 0 | 1.01× | 7→7 |
| control `byte_exact`/INTACT | 200 | 0 | **0** | 0 | 1.01× | 3→4 |

## The finding

**Attachment is necessary but not sufficient, confirmed at scale.** The lever fires hard —
telemetry on four target molecules shows `attach_return_winding_skip` at 10–16 per molecule, i.e.
the pool holds 10–16 winding-matching conformers and *every one of them is detached*. It then
promotes a better-attached conformer for 51 of the 289, and re-perception still disagrees on all
51.

That generalises the single MEDZUR case (`INTACT` and independent re-perception still differs)
from an anecdote to a measured class. It also means the `structural`/`DETACHED` bucket's 6.02
points are **not** reachable by a return-path predicate: the pool does not contain a conformer
that both holds its ligands and re-perceives correctly.

## Why it lands anyway, default OFF

The code is correct, measured, free (1.01×), and non-regressive (0 losses over 541 molecules). It
is worth keeping because it converts an unmeasured hypothesis into a standing instrument: any
future lane that widens the pool or changes the embed can flip this lever and immediately see
whether an attached conformer has become available. Deleting it would mean re-deriving the whole
apparatus to ask that question again.

## What the charter got wrong

1. **Exposure over-stated 52×.** The charter treated all 52 `DETACHED`-but-`byte_exact` molecules
   as at risk. The guard's own predicate fires on **1**. `coordination_report` ("did the donor set
   change", input vs generated) and `ligands_attached` ("did a site go empty", generated alone) are
   different tests — the 52 lose 1–3 light/ambiguous donors (H, Si, B, F), the 301 lose whole
   multi-carbon haptic groups.
2. **"Up to 5.32 pts" was an upper bound on the bucket, not on what selection can reach.** The
   reachable amount is **0**.

## Residual, shipped honestly

**POVPIA is not caught** — metal sphere intact, a hydrogen detaches and C–N reads as C=N. 7/8,
never 8/8.
