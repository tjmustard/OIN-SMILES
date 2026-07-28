# results-v0.4.13-honest — frozen authoritative table for v0.4.13

**`byte_exact` 3623 (72.46%) → 3794 (75.88%), +3.42 points, 171 molecules, 0 in a bad direction.**

## 🔴 This is a DERIVED table. No sweep was run, and that is deliberate.

Derived from `results-v0.4.8-honest` by `tools/fold_transition_sim.py --arm veto`, with the
generator **not** re-run — exactly the way `results-v0.4.8-honest` was itself derived from
`results-v0.4.6-sweep` (read its `FROZEN.md`: *"The generator was NOT re-run"*). **No generator
sweep has been run in this project since v0.4.6**; v0.4.10, v0.4.11 and v0.4.12 ran zero between
them, each riding the carry-forward licence.

v0.4.12's close-out listed "a comparable ~55 CPU-h re-sweep" as a precondition for promoting the
donor fold. That precondition was **over-specified**. What a promotion owes is a new authoritative
*table*, and here the offline route is not merely cheaper — it is **strictly better evidence**:

1. A fresh sweep re-runs a **stochastic** generator, so run-to-run variation between two
   independent runs would contaminate the +3.42 the release is reporting. The project's own
   standing trap: *never A/B by re-running a stochastic harness.*
2. It would not be like-for-like with the v0.4.8 table it is meant to be compared against, which
   is itself a re-score of a *different* generator run.

## The condition that makes it exact, which was measured and not assumed

The offline route is valid only if the promotion cannot change **what the generator returns**.
The fold is encoder-side, but that alone is insufficient: `accept_fn` (`_reencode_key_matches`)
accepts a conformer by comparing `canonical_roundtrip_key(...)` against the target key, so a fold
that moved a key would change acceptance, change the returned conformer, and make this table
describe structures the shipped code never produces.

`tools/fold_key_invariance.py`, over **every** string in the frozen corpus:

| | |
|---|---:|
| strings compared | **9669** |
| strings the fold MOVED | **1019** ← the lever is genuinely wired |
| strings whose KEY changed | **0** |
| skipped (unparseable) | **0** |

The 1019 is the load-bearing half: a lever that never toggled would *also* print 0 key changes.
Verdict **GENERATOR_NEUTRAL**. This also reproduces v0.4.11's "key untouched on 0 of 992" at ~10×
the scale.

## Denominators for the table itself

`fold_transition_sim.py --arm veto` re-encoded **393 movers from coordinates**, with
**0 excluded for drift** and **0 excluded as unavailable**. The drift control re-encodes each
molecule with the fold OFF first and requires it to reproduce the stored v0.4.8 string; passing
on all 393 independently re-confirms that v0.4.9/v0.4.10/v0.4.11/v0.4.12 left the default path
byte-identical.

⚠ An earlier run of this same arm excluded **all 393** (relative `--dataset` defaults, run from a
worktree) and still printed `+7.86 points` — the **bare fold's refuted number** — under a heading
saying "veto". The tool now refuses to report when it measures 0 of N movers.

## Contents

| file | what |
|---|---|
| `bucket_report_PASS1_authoritative.{md,json}` | **the frozen record** |
| `fold_transition_veto.json` | per-molecule transitions, all 171 named |
| `fold_key_invariance.json` | the generator-neutrality proof |
| `attach_class_audit.json` | Lane 2 — the MEDZUR / GAVSED split |
| `mirror_armA_promoted.json` / `mirror_armB_noveto.json` | the two-arm safety gate |
| `SOURCE` | pointer to the base table |

## What this table does NOT contain, and where to get it

- **Per-molecule strings for the 171.** They are identified by name in
  `fold_transition_veto.json`; the strings themselves are reproducible by re-encoding
  `results-v0.4.6-sweep/structures/<mol>_generated.xyz` with the shipped defaults.
- **Runtime.** Nothing default-path changed in generation (generator-neutral, above), so the
  v0.4.8 runtime figures stand unchanged: n = 5000, median **7.19 s**, max **759.9 s**,
  `> 30 s` **994 (19.88%)**. Recompute from the **nested** `metrics.elapsed_s` — a top-level read
  yields `None` on 5000/5000 reports — and remember it is a **SUM** over up to three separately
  SIGKILLed attempts, which is why max reads 759.9 s against a 300 s budget.

## The safety gate this table rests on

`tools/mirror_audit_donor_fold.py`, uniform n = 250, seed 7, on `cohort-v0.4.5-5k`:

| | promoted defaults | `OIN_FOLD_PARITY_VETO=0` |
|---|---:|---:|
| `REGRESSION_raw_collapsed` | **0** | **33** |
| `achiral_or_preexisting_fold` | 134 | 134 |
| `distinct_both_arms` | 116 | 83 |

**83 + 33 = 116.** The veto converted precisely the 33 collapsing molecules into separated pairs
and moved nothing else — so the zero is not bought by declining everything, which is exactly how
v0.4.12's own first veto implementation failed (it declined on 18 of 18 while all three fixture
tests passed).

🔴 **Neither `byte_exact` nor the round-trip key can see this fold's damage** —
`compare._parse_vertex_colors` folds that axis deliberately. Only the mirror audit can. Any future
release that touches `OIN_CANONICAL_DONOR_FOLD` or `OIN_FOLD_PARITY_VETO` must re-run it; the gate
arms and this table will both stay green while enantiomers collapse.
