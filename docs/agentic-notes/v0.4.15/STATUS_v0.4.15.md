# v0.4.15 STATUS — 2026-07-29, mid-release

**Base:** `main` @ v0.4.14 · **Baseline:** `results-v0.4.14-sweep`, `byte_exact` **77.16%** (n=5000)

`spec/handoffs/` is gitignored, so this file — not `spec/handoffs/v0.4.15/PROGRESS.md` — is the
record that survives.

---

## Done and committed

| | state |
|---|---|
| Baseline sweep frozen | ✅ `measurements/v0.4.14-sweep/` (bucket reports, per-molecule extract, 3 audits) |
| Charter re-derived from measured numbers | ✅ `BASELINE.md` + `BASELINE_SWEEP_CORRECTIONS_v0.4.15.md` |
| Lane 1 `OIN_ATTACH_RETURN` (+`_STRICT`) | ✅ `swimlane/v0415-attach` @ `9a1d94ca` — suite **1007 OK** |
| Lane 2 `OIN_ACCEPT_STRING_EXACT` | ✅ `swimlane/v0415-enantiomer` @ `84372966` — suite **1027 OK** |
| Combined arm | ✅ `swimlane/v0415-both` @ `63e5563e` — suite **1039 OK** |
| Populations frozen | ✅ 5 files in `measurements/v0.4.15/` (289 / 52 / 365 / 201 / 200) |
| Arm runner | ✅ `tools/run_v0415_arms.sh` |
| Lane 1 pre-flight | ✅ `measurements/v0.4.15/attach_return_preflight.json` |

All three suites green with **every lever unset**, so the default path is byte-identical. Nothing
has been pushed, per the standing instruction on this repo.

## 🔴 CORRECTION 2026-07-30: the first six arms were VOID — they measured main's tree

`generator_ab_honest.py` line 78 inserts its own `../src` at `sys.path[0]`, which **overrides
`PYTHONPATH`**. The runner set `PYTHONPATH` to the lane's worktree but invoked *main's* copy of the
tool, so all six arms imported main's `oinsmiles` — where neither lever exists — and both sides of
every A/B ran identical code. They returned a flawless `0 gains / 0 losses / output moved 0` over
1107 molecule-pairs.

Caught by `tools/selection_pool_probe.py` contradicting the telemetry, not by inspection. Fixed by
invoking `$SRC/tools/generator_ab_honest.py` plus a hard refusal that resolves
`oinsmiles.__file__` and exits unless it is under the arm's tree. Void JSONs kept at
`results-v0.4.15-arms/INVALID-wrong-tree/`. Full write-up: §9 of
`BASELINE_SWEEP_CORRECTIONS_v0.4.15.md`.

**Still valid:** every telemetry reading below, the Lane 1 pre-flight, and all three suites.
**Void:** the six arm JSONs. Arms re-launched from the correct trees 2026-07-30.

## Running

`tools/run_v0415_arms.sh lane1` and `lane2`, launched 2026-07-29 23:13, output to
`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.15-arms/`:

| arm | population | n |
|---|---|---:|
| lane1 | `pop_L1_target_site_lost` | 289 |
| lane1 | `pop_L1_control_byte_exact_detached` | 52 |
| lane1 | `pop_control_byte_exact_intact_200` | 200 |
| lane2 | `pop_L2_mirror_match` | 201 |
| lane2 | `pop_L2_target_key_equal` | 365 |
| lane2 | `pop_control_byte_exact_intact_200` | 200 |

6 single-threaded processes (`NLWP` 1 verified, BLAS capped), 74–94% CPU each on 12 cores.
Early throughput suggests **~4–8 h**. ⚠ System load average reads ~18, but that is
`npm exec ccstatusline` churn, **not** contention on the arms — per-process CPU% is the direct
evidence and a starved process would read low, not 90%.

## 🔴 PARTIAL RE-RUN RESULTS (in flight) — and the scope decision is what makes Lane 2 work

Controls **complete**, and they are clean:

| arm | n | gains | losses | runtime ratio ON/OFF |
|---|---:|---:|---:|---:|
| L1 `byte_exact`/DETACHED control | 52 | 0 | **0** | 1.008 |
| L1 `byte_exact`/INTACT control | 200 | 0 | **0** | 1.006 |
| L2 `byte_exact`/INTACT control | 200 | 0 | **0** | 1.001 |

**452 control molecules, zero regressions, and essentially zero runtime cost** — both levers are
free on the population that already passes, because those molecules accept on an early conformer
and neither predicate ever fires. The cost lands only where the pool actually fills.

Target arms:

| arm | done | gains | losses | output moved |
|---|---:|---:|---:|---:|
| L2 `key_equal` (365) | **365 ✅** | **48** | **0** | 48 |
| L2 `MIRROR_MATCH` (201) | **201 ✅** | **1** (`TAYDUV_comp_0`) | **0** | 1 |
| L1 `SITE_LOST` (289) | 250 | **0** | 0 | 46 |

### Lane 2 — MEASURED: +48 molecules = **+0.96 pts**, zero losses

The decomposition, computed against the frozen `MIRROR_MATCH` membership:

| subset | n | gains | rate |
|---|---:|---:|---:|
| **non-enantiomer** (`slot_renumber` / `rdkit_canonical`) | 164 | **47** | **28.7%** |
| enantiomer (`MIRROR_MATCH`) | 201 | 1 | **0.5%** |

**28.7% vs 0.5% — a 57× difference.** Everything the lane recovers is outside the class the
charter aimed it at. Scoped to the chartered 201 it would have returned +1 molecule.

### 🔴 And the cost is severe: 4.00× runtime on the affected population

| | OFF | ON |
|---|---:|---:|
| total elapsed over the 365 | 4294 s | **17191 s** (**4.00×**) |
| of the 365, `> 30 s` | 30 | **122** (**4.1×**) |

Extrapolated to the corpus: **+0.96 pts `byte_exact`** against **+92 molecules over 30 s**
(678 → ~770, +13.6%) and ~+3.6 CPU-h on a 38.7 CPU-h sweep. The roadmap's target is
`byte_exact` 100% **and** `max(elapsed_s) < 30 s`, so this is close to a wash between the two
halves of the goal and is **an owner decision, not a default I should pick**. It is exactly the
predicted direction (">30 s UP, materially"), at the top of the predicted magnitude.

The lane is non-regressive in accuracy (0 losses across 365 + 200 control) and the whole cost is
latency — as designed. But "free" it is not.

### 🔴 The decomposition, which is the release's real finding

**All of Lane 2's gain is in the `slot_renumber` portion of `key_equal` — none of it is in the
enantiomer class.** 9 of the first 50 `key_equal` molecules became `byte_exact` (~18%), while 0 of
45 `MIRROR_MATCH` molecules did.

That splits Lane 2 cleanly in two:

* **`slot_renumber`**: the pool *does* contain a string-exact conformer; acceptance was stopping on
  a merely-key-equal one first. A genuine selection bug, and the lever fixes it.
* **The 201 enantiomers**: the pool contains exactly **one** key-matching conformer and it is the
  mirror (telemetry §8). No acceptance predicate can fix these — **construction must**.

⚠ **Had Lane 2 been scoped to the chartered 201, it would have recovered ZERO.** The
owner-accepted widening to all 365 `key_equal` (2026-07-29) is the only reason the lane works. The
charter's framing — "the enantiomer class is the target, `slot_renumber` is a different lane" —
was measured wrong: the enantiomers are the *unreachable* part and `slot_renumber` is the
*reachable* part.

Lane 1 is moving structures (4 of 22) without yet converting one, i.e. it promotes attached
conformers that still fail re-perception — the MEDZUR shape, where attachment is necessary but not
sufficient.

**These are partials. Do not quote them as rates.**

## The earlier finding, now correctly scoped to the enantiomers only

Lane 2's lever **fires correctly and recovers nothing** on the molecules probed so far. Telemetry
(`OIN_TELEMETRY=1`) on AFADOC_comp_0 and AGAVIQ_comp_0:

```
pool.accept_incumbent_recorded     1        <- exactly ONE pool conformer carried the key
pool.accept_incumbent_returned     1        <- and it was the mirror, so it came back unchanged
adapter.string_exact_incumbent     2 / 5
```

**One key-matching conformer in the whole pool, and it is the wrong enantiomer.** No acceptance
predicate can fix that — only construction can. This is Lane 2's chartered Q4 answered *negative*
for these two, and the retain-incumbent design means the feared "183 silent wrong answers become
183 loud failures" does **not** happen: the incumbent returns and nothing regresses.

Two molecules are not a rate. **If the arms confirm it, Lane 2 is refuted as a `byte_exact` lever**
and its deliverable is a generator-capability finding for a later release. That is a legitimate
outcome — the fourth time this project would end a lane by refuting its own plan.

⚠ **Both lanes hinge on the same unmeasured quantity: does the pool contain a better conformer at
all?** Lane 1 needs an attached one, Lane 2 a correctly-handed one. A near-zero result in either is
a statement about **construction**, not about the predicate.

## Next, in order

1. **Read the six arm JSONs.** Per lane: gains, losses, `GENERATED output moved`, and the control
   arms — the controls are the only thing that can see a regression among the 3858.
2. **Decide each lever's default from its arm, separately.** A combined-only number is
   unattributable and must not be quoted as either lane's. If a lane recovers ~nothing, ship it
   default-OFF with the negative result written up rather than promoting it on mechanism alone.
3. **`both` arm** — `tools/run_v0415_arms.sh both` (Lane 1 held ON, Lane 2 varying = Lane 2's
   marginal effect). Run it **after** the single-lane arms so the box is not oversubscribed.
4. **One full 5k sweep** of whatever default is chosen (~38.7 CPU-h / ~7.7 h wall, 6 shards, BLAS
   capped, `--shard` is **1-BASED**). Verify completeness by **report count = 5000**, not exit
   status. Three full sweeps was declined, so the per-arm attribution carries a sampling caveat
   that the release note must state.
5. **`harvest_measurements.py --release v0.4.15 --from <arms> --from <sweep>`** — all `--from` in
   ONE invocation; `ALLOW` already covers `attach_*.json` and `string_exact_*.json`, but check it
   covers the arm filenames (`lane1_pop_*.json`) before trusting the total.
6. `spec/handoffs/v0.4.15/CLOSEOUT.md` verbatim: predicted-vs-actual, roadmap ×2, suite, ruff,
   version bump, tag, CHANGELOG. **Commit; do not push.** Then generate `v0.4.16/`.

## Prediction on record (CLOSEOUT §3 diffs it)

Baseline **77.16%** · `>30 s` **678** · `max(elapsed_s)` **728.8 s**.

| | predicted | note |
|---|---|---|
| `byte_exact` | +2 to +7 pts | ⚠ **now doubtful** — see the finding above; a near-zero outcome is live |
| `> 30 s` | UP, materially | both levers keep the pool filling where it used to stop |
| `max(elapsed_s)` | UP | same cause; the generator timeout is advisory |

Stating this before the arms land is the point. If it missed, the miss is the deliverable.
