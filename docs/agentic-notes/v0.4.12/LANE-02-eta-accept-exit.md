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
> G2 compares round-trip keys, and **both fold `|mc:|` by design**. An arm can therefore return
> a structurally different answer — up to and including the opposite enantiomer — while every
> gate reports "identical".

`_select_by_geometry`'s own comment has said this for two releases — *"accepting on the key
alone would hand back the wrong enantiomer while reporting success"* — and no A/B in this
project has ever measured it. This is the instrument v0.4.11 paid for: 393 molecules moved one
direction, key changed on 0 of 992, both gate arms clean, **221 enantiomers collapsed**.

G5 reports with its **denominator**, because *"0 divergent"* and *"0 divergent over 0
measured"* have printed identically in this script's history. On the first pilot it measured
**nothing** (`token_for_mol` returned `None` on all 8) and said so rather than passing.

---

## 5. 🔴 The A/B on the real population — a large speedup, and G5 FAILS

Re-run on a **quiet box** after the load-26 attempt was discarded. 12 eta molecules whose key
never matches and which exceed 30 s, `--timeout 150 --hard-cap 240 --workers 2`.

| | A default | **B lever** |
|---|---:|---:|
| `passed` | 4/12 | **4/12** |
| median | 61.8 s | **18.82 s** (−69.5%) |
| total | 824.8 s | **433.5 s** (−47.4%) |
| `> 30 s` | 11 | **5** |
| independent re-perception | 1/12 | **3/12** |
| `clash_vdw` | 0 | 1 (`PAWJED`, worst_overlap 0.7509 → 0.6073) |
| `sha256(smiles_2)` identical | — | **11 / 12** |

**The runtime thesis is confirmed**, and it is bimodal exactly as v0.4.10's speed work was —
`RIRYOJ` 119.04 s → **4.4 s** (27×, byte-identical string, passing in both arms), `PEDPOG`
42.16 → **1.93 s**, `NODLEA` 61.46 → **6.32 s**, `TEQHOM` 58.23 → **6.18 s**; while `UQUXAG`
(11.06 → 12.25) and `MOSLEL` (62.15 → 61.29) are nil.

**G2 nets positive but is not one-directional:** 3 `indep` fixes (`PAWJED`, `PEDPOG`, `NODLEA`)
against **1 regression (`KIHHUG`)**. No pass regressions, no pass fixes.

### And then the fifth arm did its job

```
G5 METAL CONFIGURATION divergent: 1 over 12 measured (of 12 molecules) ['KIHHUG_comp_0']
      KIHHUG_comp_0: A='' -> B='|mc:-|'   sha_out IDENTICAL
```

`KIHHUG_comp_0`'s two arms returned structures whose **metal-configuration descriptors differ** —
arm A's had no perceptible helicity, arm B's reads Δ/Λ-minus — and the **emitted OIN string is
byte-identical** (`5d89b34d4699566e` in both). It is also the single `indep` regression.

> **Every one of G1–G4 reports this molecule as unchanged.** G3 compares `sha256(smiles_2)`:
> identical. G4 compares pass rate: unchanged. The divergence is visible only to an arm that
> reads the descriptor `compare.py` folds by design.

This is the instrument v0.4.11 paid for, catching a real instance on its first proper run.

**Stated precisely, because overstating it would be the same sin.** An absent token means no
helicity was *perceived*, which is **not** the same claim as "the opposite enantiomer was
returned". The tool's message was corrected to print the tokens and let the reader draw the
conclusion. What is certain: the arms returned **structurally different** answers on a molecule
that every existing gate calls identical, and that is disqualifying on its own.

**Verdict: `OIN_ETA_ACCEPT_EXIT` stays default-OFF.** The charter made any Δ/Λ divergence a
*blocking* finding, and one appeared at n=12. The speedup is real and large; it is not free, and
the thing it costs is invisible to four of the five gates.

## 6. What is deliberately NOT measured here

**Nothing, in the end — but the first attempt was thrown away and that was right.** The A/B
was initially run at load 26 with two mirror audits going, and the standing trap is *never
interleave timing runs with gate runs*. It was **stopped rather than banked** and re-run on a
quiet box. The two runs disagree by enough to have wrecked the conclusion: `UQUXAG_comp_0`
reads **17.93 s** under load and **11.06 s** clean — a 62% inflation, comparable to the entire
effect being measured. Under load the per-molecule hard cap also fires spuriously and
manufactures timeout-shaped "regressions" (v0.4.4 read 11 of those as correctness deltas; all
11 were timeouts).

Also out of scope, and handed forward in writing: the **GAVSED class** —
`_select_by_geometry`'s fallback ranking is not attachment-aware, so the check guards
*acceptance*, not *return*. Closing it changes arm A's behaviour and needs its own gate.

---

## 7. Files

`src/oinsmiles/generation/metallogen_adapter.py` (`_geometry_classifies`,
`_eta_accept_exit_ok`, `accept_fn` wiring) · `src/oinsmiles/oin/levers.py` ·
`tools/ab_accept_scored.py` (`--lever`, `--extra-env`, G5) ·
`tests/unit/test_eta_accept_exit.py` (16 tests).
