# v0.4.13 — the donor fold ships, and four instruments were caught printing nothing

> **`byte_exact` 72.46% → 75.88%, +3.42 points, 171 molecules, 0 in a bad direction.**
> The first time this project's headline has moved **up**. The only other time it moved at all was
> v0.4.8, which took it *down* 10.34 points on purpose, to make it honest.

Lane detail: [`LANE-02-attach-classes.md`](LANE-02-attach-classes.md) ·
table `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.13-honest/`

---

## 1. What shipped

`OIN_CANONICAL_DONOR_FOLD` and `OIN_FOLD_PARITY_VETO` go default-ON **together**. They are one
promotion with two names: the fold is the change that pays, the veto is the condition under which
it is safe, and the pair is the smallest unit that is both.

v0.4.11 built the fold, measured +7.86 points across 393 molecules, and refuted it — it collapses
enantiomers in **221 of those same 393 gains**. v0.4.12 built the reflection-parity veto and
measured 171 survivors. This release runs the gate on two populations and ships it.

| | cat/ (v0.4.12's population) | cat+photo (broader) |
|---|---:|---:|
| `REGRESSION_raw_collapsed`, veto OFF | **19** | **33** |
| `REGRESSION_raw_collapsed`, promoted | **0** | **0** |
| `achiral_or_preexisting_fold` | 157 → 157 | 134 → 134 |
| `distinct_both_arms` | 73 → 92 | 83 → 116 |
| accounting | **73 + 19 = 92** ✓ | **83 + 33 = 116** ✓ |

The accounting line is the point. A veto that works and a veto that declines everything both
produce a separation count of zero collapses — v0.4.12's own first implementation declined on
**18 of 18** while all three fixture tests passed. Here every collapsing molecule moved into
`distinct_both_arms` and **nothing else moved at all**, so the zero is demonstrably bought by
separating rather than by abstaining.

The cat/ arm reproduces v0.4.12's published table (`157 / 73 / 1`) line for line. The broader draw
reads 13.2% against cat-only's 7.6% — a second population, not a contradiction, and it says the
collapse rate is not an artifact of the corpus half v0.4.11 happened to sample.

## 2. 🔴 The finding that generalizes: the 55 CPU-h sweep was the wrong instrument

v0.4.12's close-out made a full re-sweep a precondition for this promotion. **That precondition
was over-specified, and following it would have produced weaker evidence at 9 hours' cost.**

Three facts, in the order they matter:

1. **`results-v0.4.8-honest` — the 72.46% baseline — is not a sweep.** It is an offline re-score of
   `results-v0.4.6-sweep`. Its own `FROZEN.md`: *"The generator was NOT re-run."*
2. **No generator sweep has been run in this project since v0.4.6.** v0.4.10, v0.4.11 and v0.4.12
   ran zero between them.
3. **The generator is stochastic.** A fresh run differs from the frozen corpus for reasons that
   have nothing to do with the fold, and the standing trap is *never A/B by re-running a
   stochastic harness*. A sweep would have contaminated a +3.42-point signal with run-to-run
   variation, and would not even have been like-for-like with the table it was compared to.

What a promotion owes is a new authoritative **table**, not a new generator run.

### The condition that makes the offline route exact — measured, not assumed

The fold is encoder-side, but that alone is not sufficient. `accept_fn`
(`_reencode_key_matches`) accepts a conformer by comparing `canonical_roundtrip_key(...)` to the
target key, so a fold that moved a key would change acceptance, change the returned conformer, and
make an offline table describe structures the shipped code never produces.

`tools/fold_key_invariance.py`, over **every** string in the corpus:

| | |
|---|---:|
| strings compared | **9669** |
| strings the fold MOVED | **1019** |
| strings whose KEY changed | **0** |
| skipped | **0** |

The 1019 is the load-bearing half: **a lever that never toggled would also print 0 key changes.**
Verdict `GENERATOR_NEUTRAL`, and it reproduces v0.4.11's "0 of 992" at ~10× the scale.

## 3. Three instruments, caught printing plausible nothing

This release's own rule — *before quoting an instrument, ask what a BROKEN version would print* —
fired three times in one session, twice on instruments this project already trusted.

**(a) The mirror audit died on a deleted dataset.** A sibling session checked out the data-only
`Kulik_TMC_Dataset` branch in the main checkout; switching back to `main` **deleted all 26,232
files** under `cat/`+`photo/`, because they are tracked there and untracked on `main`. Every cohort
symlink dangled. Recovery is `git restore --source=Kulik_TMC_Dataset --worktree -- <paths>`, and
the standing check is `find <cohort> -xtype l | wc -l` **= 0** before trusting any measurement.

**(b) `fold_transition_sim.py` printed the REFUTED number.** Run from a worktree, `--dataset`
defaults to a **relative** path that exists only in the main checkout, so `_find_input_xyz`
returned `None` for all 393 movers, every one landed in `unavailable`, and the run printed:

```
⚠ EXCLUDED 393 molecules: structure not on disk
byte_exact 3623 -> 4016  (72.46% -> 80.32%, +7.86 points)
```

and exited 0. **That +7.86 is the bare fold's number — the one v0.4.11 refuted for collapsing
enantiomers — printed under a heading that says "veto".** The mechanism is the tool's own correct
choice to let excluded movers keep their fold-arm bucket: pessimistic for a *partial* exclusion,
catastrophic for a *total* one, where the veto arm degenerates into the fold arm.

This is worse than the usual shape. The broken version did not print the same thing as the working
one; it printed something *more attractive* and wrong. The tool now refuses when it measures 0 of
N movers, and prints `veto arm measured M/N` when it succeeds.

**(c) `run_sweep.sh` would have run on the wrong interpreter.** It still carried the sibling-glob
fallback that `gate_v047.sh` was hardened against in v0.4.9 — live in the one place it costs 55
CPU-h. From the v0.4.13 worktree the glob selected `EtaCatalysis/.venv` (no rdkit at all); the next
candidate is `EtaTMCSMILES/.venv` at rdkit **2025.09.2** against the pinned **2025.9.3** — the loud
failure and the silent one, in that order. Now resolves via `--git-common-dir` and refuses on drift.

## 4. The classes that were n = 1 are 280 and 99

Full detail in [`LANE-02-attach-classes.md`](LANE-02-attach-classes.md). Over the **767 genuine
failures** (excluding `key_equal`, which is benign canonicalization and not a disagreement):

| | n = 1 → |
|---|---:|
| **GAVSED shape** (`DETACHED` — returned with ligands off the metal) | **280** |
| **MEDZUR shape** (`INTACT` — attachment fine, re-perception still disagrees) | **99** |

Control: `byte_exact` reads **1.32%** `DETACHED` against **24.11%** on the failing side — **18.2×
enrichment**, `UNKNOWN` = 0, `#DONE 5000`.

**🔴 It re-mechanises a release four rungs up.** `structural` — 417 molecules, **8.34 points**, the
second-largest block in the gap, scheduled v0.4.17 as *"bounded by what the generator can
assemble"* — is **266/417 = 63.8% `DETACHED`**. Those did not fail from missing capability. The
generator assembled a structure and **returned** it with ligands off the metal, because
`_select_by_geometry`'s fallback ranking is not attachment-aware. The guard already exists for
*acceptance* (`OIN_ATTACH_CHECK`); applying it to *return* is one site.

## 5. 🔴 A fourth blind instrument: ARM 1 cannot see this release at all

`tools/gate_v047.sh arm1` returns **PASS — byte-identical to golden, `#DONE 62`**, with
`MANIFEST_SHA256=6c2417e1…` — **the same hash before and after the promotion.**

The tempting reading is "the gate validates the promotion." It does not. Checking the 62 fixture
names against the 393 fold-movers:

```
fixtures in ARM 1:               62
fixtures that are fold-movers:    0
```

**ARM 1 has zero coverage of the population this release moves.** Its PASS means "no regression on
the fixtures", which is worth having, but it is *not* evidence that the fold works, and a session
that read it as such would be repeating this release's central mistake with a different tool. The
evidence for the promotion is the mirror audit and the transition simulation, nothing else.

The same check on the round-trip arm found the opposite, which is why it was worth running:

| golden | rows | fold-movers |
|---|---:|---:|
| `gate_v047_arm2_golden.tsv` | 100 | **0** |
| `gate_v049_arm2_golden.tsv` | 325 | **11** |

So ARM 1 and the v0.4.7 ARM 2 golden need **no** re-freeze — not as an economy, but because the
promotion provably does not touch them — and exactly 11 rows of the v0.4.9 golden do. Those 11 were
re-run individually and patched, rather than bulk-accepting a whole regenerated manifest, so the
diff stays reviewable.

**Standing consequence:** the byte-identity gate is structurally blind to canonicality levers whose
movers are corpus molecules rather than curated fixtures. Any future release promoting one must
state its gate's *coverage of the moved population*, not just its verdict.

## 6. Lane 1 — the prefilter defect is real, and overriding it recovers nothing here

`OIN_PREFILTER_ADVISORY` is built, registered default-OFF, and instrumented. The wiring gate on
`AROHIA_comp_0` — where the answer is known independently — reads, **scored honestly**:

| | lever OFF | lever ON |
|---|---|---|
| `prefilter_veto_overridden` / `_confirmed` | — | **2** / 1 |
| how acceptance ended | `early_exit_`**`miss`** | `early_exit_`**`hit`** |
| what was returned | a **previously-rejected** conformer, via the geometry fallback | a conformer the strict test **accepted** |
| round-trips (honest) | **True** | **True** |
| elapsed | 22.73 s | **2.52 s** |

The defect is **confirmed**: the cheap prefilter rejects conformers the strict test takes. But
**`recovered = 0` because both arms already pass**, not because both fail — the accuracy delta on
this molecule is nil. What changes is *how* the pass is obtained: lever-off, `accept_fn` accepts
nothing and `_select_by_geometry` returns a conformer acceptance had already rejected — the same
unguarded fallback Lane 2 measures at **280 molecules corpus-wide**. Lever-on, a conformer is
accepted on its merits, 9× faster.

Three cautions, all of which change how the numbers read:

1. **The denominators are not pool sizes.** The first override *accepts*, which stops the pool
   filling, so the lever-ON arm evaluates far fewer conformers than the pool would hold. AROHIA's
   documented 0/48-vs-16/48 was measured with the pool **forced full**; "3 vetoes" here is the same
   defect observed until it stopped mattering, not a smaller one.
2. **The negative latency delta (−20 s) is not a speed claim** — it is early exit, and it was
   measured on a loaded machine besides.
3. **🔴 The first version of this tool scored with `get_oin_string(res.mol, coords)` — the
   circular predicate v0.4.8 replaced — and it INVERTED the answer.** Identical telemetry, opposite
   verdict: circular reads `passed=False` in **both** arms, honest reads `passed=True` in **both**.
   Anything written from that run would have said "the molecule fails either way"; the truth is "it
   passes either way, by different routes". A live single-molecule instance of the 8 false
   negatives v0.4.8 measured, reproduced by accident while building something else — and a
   reminder that the circular predicate still lurks in the older A/B tools this one was copied from.

**Corpus prevalence is NOT measured.** n = 1 is exactly what the charter forbids quoting, so this
lane is honestly incomplete and is handed to v0.4.14 with its instrument built and gated.

## 7. Predicted vs actual

| | predicted | actual | |
|---|---|---|---|
| `byte_exact` | +3.2 to +3.6 pts | **+3.42 (72.46% → 75.88%)** | ✓ 171 molecules |
| mirror audit, promoted | 0 | **0 on both populations** | ✓ |
| mirror audit, veto off | 19 | **19 (cat/) and 33 (cat+photo)** | ✓ + a second draw |
| `facmer_divergent` | does not rise | **16 → 16** | ✓ |
| bad-direction moves | 0 | **0** | ✓ |
| MEDZUR / GAVSED | "could be 2 or 200" | **99 / 280** | — the sketch's deletion clause does not fire |
| the sweep | ~55 CPU-h | **not run, and should not be** | ✗ **the plan was wrong** |

The last row is the release's real deliverable. Four of this project's releases have ended by
refuting their own plan; this is the fifth, and the first to refute a *methodological* precondition
rather than a technical one.

## 8. What v0.4.14 inherits

1. **Lane 1 is built and unmeasured.** `OIN_PREFILTER_ADVISORY` (default OFF) makes the cheap
   acceptance prefilter advisory, with `overridden` / `confirmed` / `cheap_pass` telemetry and
   `tools/prefilter_prevalence.py`. The AROHIA two-point fixture is the wiring gate. **It needs a
   quiet machine** — its deliverable includes a latency cost.
2. **The GAVSED guard.** 266 molecules of `structural`, up to 5.32 pts, one site. Needs its own
   gate: closing it changes arm A's behaviour.
3. **The 48 `byte_exact` molecules that read `DETACHED`.** Two readings survive and neither was
   tested. Named in Lane 2 §6 rather than guessed at.
4. **The carry-forward licence is spent and re-earned.** The drift control re-encoded all 393
   movers with the fold OFF and required them to reproduce the frozen v0.4.8 strings. They do —
   so v0.4.9 through v0.4.12's byte-identity claims are now checked against the corpus rather than
   inherited, for the mover population at least.
