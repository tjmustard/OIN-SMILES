# Lane 2 — `OIN_ETA_ACCEPT_EXIT`: the winding criterion, moved to where it can pay

**Status: BUILT, default-OFF. The chartered lane was overridden, and the chartered cohort was
stale.** Both are results, not detours.

---

## 1. The charter asked for an A/B that could not have told us anything

The v0.4.12 sketch chartered this lane as:

> *the `OIN_ETA_EARLY_EXIT` corpus A/B its own promotion gate demands, and that has never been
> run*

Reading the code before running it:

```
metallogen_adapter.py:1715   _eta_early_targets  ->  inside _select_by_geometry_impl
metallogen_adapter.py:2139   accept_fn           ->  passed INTO generate_3d_structures
```

`_select_by_geometry_impl` runs **after** `generate_3d_structures` has already filled the
entire pool. `accept_fn` is the only site consulted per conformer **during** filling. So
`OIN_ETA_EARLY_EXIT` cannot reduce embed count by construction — and its own in-code A/B
already recorded exactly that: *Ferrocene, lever off → 32 attempts; lever on → 32 attempts. It
fires, and the attempt count does not move.*

> **Its promotion gate is VOID, not unrun.** Spending a corpus A/B to re-measure a documented
> structural null buys nothing. The lever is kept only as the marker of where the boundary is,
> and its `_HELD_OFF` entry now says so.

The criterion was relocated to `accept_fn` as `OIN_ETA_ACCEPT_EXIT`. Eta is **53.1% of the
whole `> 30 s` tail** (528 of 994), which is why this is where the runtime goal lives.

---

## 2. It is a conjunction, and that is not defensive padding

Accepting on winding **alone** would stop the pool before `_select_by_geometry`'s clash-first
ranking ever ran — structurally the *same* defect `OIN_ACCEPT_SCORED` has, which cost **26
independent-re-perception regressions on a 100-molecule population**. So a conformer must also:

1. classify as the requested coordination geometry (`classify_and_fit`, hoisted to
   `_geometry_classifies`), and
2. still have every claimed coordination site populated (`conformer_ligands_attached`).

The attachment check is **unconditional inside this branch** rather than gated on
`OIN_ATTACH_CHECK`. v0.4.7's one unambiguous finding was *never run a scored-acceptance lever
without it*; here it is part of what the predicate **means**, and a test asserts that no lever
read stands between the branch and the check.

Anything unevaluable returns `False` — the pool then fills exactly as it does today. Precedence
is unchanged: an exact key match still wins, so a molecule that round-trips today is untouched.

---

## 3. 🔴 The chartered cohort no longer describes the acceptance gap

The obvious pilot population is `docs/agentic-notes/v0.4.6/eta_accept_gap_cohort.md`'s
`CHEAP_ONLY` class — molecules where the strict test *never* matches, so the pool always fills
(`GAVSED` 319.9 s, `HEJXIF` 147.6 s, `WIWRIE` 77.6 s, `NOMMOU` 52.5 s). An A/B over its 8
`GAP`/`CHEAP_ONLY` molecules came back **completely flat**:

| | A default | B lever |
|---|---:|---:|
| median | 40.64 s | 40.86 s |
| total | 418.1 s | 421.6 s |
| `> 30 s` | 5 | 5 |
| `sha256(smiles_2)` identical | — | **8 / 8** |

A flat A/B is only a result if you can name the mechanism. Telemetry named it, and it is *not*
that the lever is broken:

```
GAVSED_comp_0, lever off:  {"adapter.early_exit_hit": 1}
GAVSED_comp_0, lever on:   {"adapter.early_exit_hit": 1}
```

`adapter.eta_accept_*` never fires **because the key path already succeeds** — the fallback is
only reached when the key check fails. All 8 carry live eta targets, so the branch is active;
they simply are not in the gap any more. **The v0.4.6 cohort predates v0.4.8–v0.4.11 and its
classification is stale.**

### The real target population, re-derived from the frozen honest sweep

```bash
# eta molecules = smiles_1 matching \{\d+[<>]\}, bucketed with roundtrip_bucket_report.classify
```

| bucket | n |
|---|---:|
| `byte_exact` | 572 |
| `structural` | 254 |
| `key_equal` | 169 |
| `hard_fail` | 138 |
| `facmer_divergent` | 13 |
| **total eta** | **1146** |
| **key does NOT match** | **405** |
| **…of those, `> 30 s`** | **378** |

So the lever's population is **378 eta molecules that today burn the budget and fail**, not the
handful of fast round-trippers the stale cohort pointed at. Its effect is to convert pool-fill
timeouts into decided outcomes — which is a *different claim* from "speed up passing
molecules", and it is gated accordingly.

---

## 4. The fifth gate arm, and why four were not enough

`tools/ab_accept_scored.py` was hardcoded to `OIN_ACCEPT_SCORED`; everything else in it is
lever-agnostic, so it gained `--lever` / `--extra-env` rather than being forked. It also gained
**G5**:

> **G1–G4 are structurally blind to metal configuration.** G3 compares `sha256(smiles_2)` and
> G2 compares round-trip keys, and **both fold `|mc:|` by design**. An arm can therefore hand
> back the **opposite enantiomer** and every gate reports "identical".

`_select_by_geometry`'s own comment has said this for two releases — *"accepting on the key
alone would hand back the wrong enantiomer while reporting success"* — and no A/B in this
project has ever measured it. This is the instrument v0.4.11 paid for: 393 molecules moved one
direction, key changed on 0 of 992, both gate arms clean, **221 enantiomers collapsed**.

G5 reports with its **denominator**, because *"0 divergent"* and *"0 divergent over 0
measured"* have printed identically in this script's history. On the first pilot it measured
**nothing** (`token_for_mol` returned `None` on all 8) and said so rather than passing.

---

## 5. What is deliberately NOT measured here

**The timing arms.** The box was at load 26 with two mirror audits running, and this project's
standing trap is *never interleave timing runs with gate runs* — a `pkill -f` incident once
drove load to 35 and corrupted every in-flight timing. Under load the per-molecule hard cap
also fires spuriously, manufacturing timeout-shaped "regressions" (v0.4.4 read 11 of those as
correctness deltas; all 11 were timeouts). The A/B was **stopped rather than banked**.

Also out of scope, and handed forward in writing: the **GAVSED class** —
`_select_by_geometry`'s fallback ranking is not attachment-aware, so the check guards
*acceptance*, not *return*. Closing it changes arm A's behaviour and needs its own gate.

---

## 6. Files

`src/oinsmiles/generation/metallogen_adapter.py` (`_geometry_classifies`,
`_eta_accept_exit_ok`, `accept_fn` wiring) · `src/oinsmiles/oin/levers.py` ·
`tools/ab_accept_scored.py` (`--lever`, `--extra-env`, G5) ·
`tests/unit/test_eta_accept_exit.py` (16 tests).
