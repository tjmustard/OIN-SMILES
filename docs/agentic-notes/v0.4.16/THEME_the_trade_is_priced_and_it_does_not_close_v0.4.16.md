# v0.4.16 — the trade is priced, and it does not close

**FLAT at 77.16%.** Both lanes delivered what they were chartered to deliver; neither moved the
headline, and in both cases that is the result rather than a shortfall.

| | predicted | actual |
|---|---|---|
| `byte_exact` | UP 0.0–0.96 | **0.00** — the bar is unreachable, so the lever stays default-OFF |
| `> 30 s` | unchanged unless promoted | unchanged |

The charter said the headline's *ceiling* was known before the release started and that "the open
question is the runtime price, not the gain". That framing was right, and the answer to the open
question is **the price does not come down**.

---

## Lane 1 — `OIN_ACCEPT_STRING_EXACT`: HOLD, because there is no knee

v0.4.15 left +48 molecules (+0.96 pts, 0 losses) sitting behind a **4.00×** runtime cost, on the
theory that the cost was a reclaimable tail: the lever declines to stop the pool, and the 317
molecules that never gain consume **93.9%** of its bill. Lane 1 built the bound and priced it.

**The two halves of the bar are mutually unreachable:**

| | needs | and then gives |
|---|---|---|
| ≥ 75% of the 48 recovered | bound **≥ 12** | `> 30 s` = **104** (limit 52) |
| `> 30 s` ≤ 52 | bound **≤ 3** | recovered **19** of 48 |

**The shape is the finding, not the threshold.** To keep 79% of the gain you pay 89% of the
runtime penalty. The frontier is close to linear, so bounding moves *along* the v0.4.15 trade
rather than improving it. The charter's mechanism hypothesis is refuted: the cost is not a wasted
tail — each post-incumbent conformer is a full embed **plus** a full `XYZToSMILES` re-encode, so
the *early* extra conformers are the expensive ones. (Charter candidate 1, "check the pool is not
filling past the first string-exact conformer", was already true: `early_hit` breaks both loops.)

Full curve, validation and the two points the owner may wish to override:
`LANE-price-the-string-exact-trade.md`.

## Lane 2 — the 3.74 "unexplained" points: **82% is PERCEPTION**

All 187 classified, 0 unaccounted, every class read. The block the roadmap carried for three
releases as *"attachment fine, re-perception still disagrees, nobody knows why"* is, over the 172
`structural` INTACT+BOUNDARY:

**PERCEPTION 141 (82.0%) · CONSTRUCTION 25 (14.5%) · STEREO 6 (3.5%)**

This re-points the ladder: v0.4.17 keeps `DETACHED` (301 / 6.02) as construction and **declines**
the 172, and the 141 perception molecules (2.82 pts) need a lane that does not exist yet.
`facmer_divergent` is 100% `ARRANGEMENT_ONLY` — the one block on the board whose name was a
measurement rather than a hypothesis. Details: `LANE-characterise-the-unmeasured.md`.

---

## 🔴 What this release is actually worth: three instruments that would have lied

Neither lane moved a point. Both produced a number that changes what the next two releases do —
and in **three separate places** the first version of an instrument produced a plausible, precise,
self-consistent, wrong answer. Recorded in full in
`METHOD_one_run_beats_a_parameter_sweep_v0.4.16.md`.

1. **A normalizer wrong 57 times in 109.** String-level heavy-atom normalization vs an RDKit
   canonical comparison: `SKELETON` 74/172 (43%) against 4/172 (2.3%). **A factor of 18.** The
   broken version would have *confirmed* the roadmap's existing assumption that `structural` is
   construction work — which is exactly why it would not have been questioned. Caught because a
   read example (`BEDLII_comp_0`) disagreed with the instrument. Neither alone sufficed: the
   eyeball read was too shallow on a 300-character macrocycle, and the instrument was confident.

2. **A derived runtime curve that understated every bounded row.** The cost model used the in-loop
   telemetry stamp and ignored the post-loop tail — selection, return, write-out. The error grew
   with the bound, i.e. it was worst exactly where the decision sits, and it biased toward
   *promoting*. Caught by an arithmetic identity the fix makes checkable: at a bound past a
   molecule's last evaluation, derived cost must equal measured `elapsed_s`.

3. **`pgrep -f` matching its own `bash -c` body** — a standing trap in this project's own notes —
   reported three live workers for an arm that had been dead for minutes, and cost a
   confirmation arm that was silently SIGTERMed at 30 of 48 without writing its JSON.

## The method result worth carrying forward

**The knee curve is arithmetic, not a parameter sweep.** `incumbent_hit` is returned whatever the
pool does afterwards, so bounding at *N* changes the answer only for molecules whose hit lies
beyond *N*. Recording each molecule's `min_bound` makes both curves derivable from **one** run —
the charter had budgeted 1–2 h *per point*.

That is only legitimate because it was validated rather than assumed, at points the derivation
could not fake:

| check | derived | independently measured |
|---|---:|---:|
| ceiling | **48** | **48** (v0.4.15 gains) — exact |
| bound-0 runtime | 4534 s | 4294 s (+5.6%) |
| unbounded runtime | 16418 s | 17191 s (−4.5%) |
| live bound-0 gate | — | **0 gains, 0 losses, 0 output moved** over 40 |

Bound 0 does not merely reach lever-OFF's *verdict*; it produces **byte-identical output**. A
broken bound and a working one agree at large *N* and differ at 0 — that is the check that
separates them.

## Durability

- `tools/freeze_sweep_extract.py` — a whole 5000-molecule sweep in **0.26 MB** (vs 268 MB for the
  directory), verified to reproduce the authoritative bucket report exactly. It carries no
  geometries and says so where it prints.
- `harvest_measurements.py`'s `ALLOW` was one pattern too narrow **for the fourth consecutive
  release**; `knee_*.json` and `pop_*.txt` would have been dropped silently. Verified with
  `fnmatch` over the real filenames, not by reading the list.

## No sweep is owed

The shipped default did not change, so the mandatory 5k sweep does not apply — the same reason
v0.4.15 skipped its own. The baseline of record remains `results-v0.4.14-sweep`, 77.16%.
