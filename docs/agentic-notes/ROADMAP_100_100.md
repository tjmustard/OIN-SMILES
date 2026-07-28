# Roadmap — `byte_exact = 100%` and `max(elapsed_s) < 30 s`

**Git-durable counterpart of `spec/handoffs/roadmap-100-100/README.md`**, which is gitignored along
with the rest of `spec/handoffs/`. The handoff tree carries the executable charters; this file
carries the ladder, so it survives a clone. **Where the two disagree, re-measure — do not pick one.**

**Authored:** 2026-07-26, against `main` @ `0c21729c` (197 commits ahead of `origin/main`, 0 pushed).
⚠ `main` moved **three times** while this was being written. Re-read the tip.
**Corpus of record:** `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.6-sweep`, N = 5000, seed-42,
`--mol-timeout 300`, PASS 1 + PASS 2 complete.

---

## Why a ladder instead of one push

Both goals are one goal, and neither is currently measurable.

1. **The headline is ~5.7 points over-stated.** The harness scores with
   `get_oin_string(gen_result.mol, coords)` — the generator's own bond graph — so it asserts bonds
   the geometry does not support. FP **61/633 = 9.6%** overall, **48/171 = 28.1% haptic**; FN 8/302 =
   2.6% the other way. `FIYHUT_comp_0` ships both Cp rings 0.85 Å off the iron (10 bonded carbons →
   0) and scores a pass. *(`docs/agentic-notes/v0.4.6/METRIC_FALSE_POSITIVES.md`)*
2. **78.8% of failures never test the notation.** Of 340 failures, 240 are generator timeouts and 28
   produced nothing. The notation-attributable gap is **57/5000 ≈ 1.1%** — where the notation is
   exercised it is ~98.9% correct. *(`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md` §2)*
3. **Speed alone does not buy accuracy honestly.** `OIN_ACCEPT_SCORED` "recovered" 90 of those 340;
   `report["coordination"]` showed **60 (66.7%) had ligands off the metal**, 21 boundary, 9 genuine.
   +90 → 95.0% is really **+9 → 93.4%**. *(`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md` §3)*
4. **The budget is not a budget.** `OIN3DGenerator(timeout=)` is advisory: 60 s requested,
   60.0–172.8 s spent. `p100 < 30 s` is not currently expressible.

---

## Definitions of done

| goal | metric | today | target |
|---|---|---|---|
| **A — accuracy** | `byte_exact`, scored by **independent** re-perception of the generated XYZ | **72.46%** *(honest, v0.4.8)* | **100%** |
| **B — speed** | per-molecule wall-clock against an **enforced** budget | **994/5000 = 19.88%** over 30 s; median 7.19 s | **< 30 s, p100** |

⚠ Two live traps in one field. `metrics.elapsed_s` is **nested** — read from the top level it
silently yields `0` — **and it is a SUM** over up to three separately SIGKILLed harness attempts.
The old "max 759.9 s against a 300 s budget" headline was the second trap: all 4658 single-attempt
rows finish within **0.2 s** of their cap. See `v0.4.9/ELAPSED_S_IS_A_SUM_v0.4.9.md`.

## The gap — `100 − 72.46 = 27.54` points (honest, re-derived v0.4.9)

⚠ **The `release` column below is post-`LADDER DECISION 2026-07-27` and post-v0.4.11.** It replaces
an earlier copy that still pointed `slot_renumber` at v0.4.14 and `encode_fail` at v0.4.11.

| block | n | pts | nature | release |
|---|---:|---:|---|---|
| `key_equal` → `slot_renumber` | 496 | 9.92 | canonicality, encoder-side | 🔴 **v0.4.11 attempted — REFUTED.** 377 (7.54) are reachable by a within-fragment fold, but that fold **collapses enantiomers in 221 of its 393 gains**. Needs a reflection-parity filter. 90 mol / **1.80 pts** of the residue re-filed to **v0.4.14** (frozen resonance form) |
| **`structural`** | **417** | **8.34** | generator capability (re-labelled) | v0.4.17 |
| `hard_fail` | 319 | 6.38 | compute (mostly) | v0.4.12 – v0.4.13 |
| `key_equal` → `rdkit_canonical` | 114 | 2.28 | canonicality, encoder-side | v0.4.14 (**+1.80 from above ⇒ ~4.08**) |
| `facmer_divergent` | 16 | 0.32 | wrong isomer | v0.4.15 |
| `encode_fail` | 15 | 0.30 | encoder coverage | v0.4.18 (opportunistic) |
| **sum** | **1377** | **27.54** ✓ | | |

### 🔴 The honest number REORDERED the ladder

`structural` was **9 molecules / 0.18 points** under the scored metric and was parked at v0.4.16 as
*"knowledge, not points"*. Honestly scored it is **417 molecules / 8.34 points — the second-largest
block in the gap**, ahead of `hard_fail`. It is not a rounding change: 350 of them came directly
from `byte_exact`, molecules whose generated structure re-perceives with different coordination
than the OIN claims.

**v0.4.16 is therefore mis-sized and mis-ordered in the ladder below, and the "wrong isomer" label
is doing too much work** — a molecule that loses a haptic contact is not an isomer problem. A
release that owns 8.34 points cannot sit behind two canonicality releases worth 12.20 combined
without someone having decided that deliberately. **Re-sequence before planning v0.4.11.**

#### ⚠ Still not re-sequenced at v0.4.10's close-out — and now one release more urgent

v0.4.10 was byte-identical by construction, so it moved **no** point of the gap and the table above
is unchanged. That means the re-sequencing decision v0.4.9 asked for was simply not taken, and the
ladder still spends its **next** release on its **smallest** block:

| block | pts | currently scheduled |
|---|---:|---|
| `key_equal` → `slot_renumber` | **9.92** | v0.4.14 |
| `structural` | **8.34** | v0.4.16 |
| `hard_fail` | 6.38 | v0.4.9–13 |
| `rdkit_canonical` | 2.28 | v0.4.15 |
| `facmer_divergent` | 0.32 | v0.4.16 |
| **`encode_fail`** | **0.30** | **v0.4.11 — next** |

**18.26 of the 27.54 points sit at v0.4.14 and v0.4.16**, behind a release worth 0.30.

There is a real argument for the current order — v0.4.12 (honest acceptance) and v0.4.13 (harness
false negatives) are *measurement* prerequisites, and moving a 9.92-point release ahead of them
risks attributing a harness artefact to a canonicality fix. That argument covers v0.4.12/13. **It
does not cover v0.4.11**, whose 0.30 points depend on nothing and unblock nothing.

**Recommendation, for the human to accept or reject** — this reorders eight planned releases and
the target definition is the user's call, so it is written down rather than applied:

1. Keep v0.4.12 and v0.4.13 where they are; they are prerequisites, not point-scorers.
2. **Re-point v0.4.11 at `slot_renumber` (9.92 pts)** — encoder-side canonicality, the largest single
   block, and the area where v0.4.5 already built six working levers, so it has the shortest path
   from diagnosis to change.
3. Demote **encode floor R3** to an opportunistic slot. It is genuinely cheap and this session hit
   its cost first-hand — ARM 1 takes ~15 minutes almost entirely because several fixtures route
   through `_resonance_candidates_isolated` at a **120 CPU-s budget per fork** — but 0.30 points is
   a rounding error against a 27.54-point gap, and "it also speeds up our own gate" is a developer
   convenience, not a roadmap argument.
4. Re-label `structural` before scheduling it. 417 molecules whose generated structure re-perceives
   with **different coordination than the OIN claims** is a generator-capability problem sitting
   under an "isomer" heading; it likely belongs with v0.4.18, not beside `facmer_divergent`.

Runtime, from the primary reports: overall n = 5000, median **7.19 s**, `> 30 s` **994 (19.88%)**;
**eta** n = 1146, median 24.08 s, `> 30 s` 528 (46.1%) — **53.1% of the whole tail**; non-eta
n = 3854, median 6.18 s, `> 30 s` 466 (12.1%).

### What the tail actually costs, and what bounding it costs (v0.4.9)

| band | n | CPU-h | % of sweep | honest passes |
|---|---:|---:|---:|---:|
| `< 30 s` | 4006 | 8.67 | 15.8% | 3380 (84.4%) |
| `30–60 s` | 292 | 3.45 | 6.3% | 136 (46.6%) |
| `60–300 s` | 411 | 15.36 | 28.0% | 112 (27.3%) |
| **`≥ 300 s`** | **291** | **27.31** | **49.8%** | **3 (1.0%)** |
| total | 5000 | **54.8** | | 3631\* |

**Half the sweep's compute buys three molecules.** And **93.1% of all honest passes already finish
under 30 s**, so the two goals are far more separable than thesis 3 below assumed:

| per-molecule cap | CPU-h recovered | passes lost | `byte_exact` cost |
|---:|---:|---:|---:|
| 300 s | 3.06 | 3 | 0.06 pts |
| 120 s | 21.33 | 52 | 1.04 pts |
| 60 s | 30.97 | 115 | 2.30 pts |
| **30 s** | **37.84** | **251** | **5.02 pts** |

A 5k sweep under a 30 s bound costs **~17 CPU-h instead of 55**, which makes live sweeping an
affordable instrument again rather than something v0.4.8 had to route around with an offline
re-score.

\* 3631 counts string equality; the frozen report says **3623**. The eight-molecule difference is
the **atom-count gate** — molecules that encode byte-identically to their input while having lost
atoms. A string comparison cannot see them. Use `tools/roundtrip_bucket_report.py`, never an ad-hoc
`honest_class.endswith("->byte")`.

---

## Ordering thesis

1. **Instrument before engine.** The headline goes **down once, on purpose**, inside its own release
   boundary — otherwise it is indistinguishable from a regression, a confound this project has been
   caught by twice.
2. **Bound before optimize.** A p100 target against an unbounded budget is not a target.
3. **Compute is the accuracy work** for the failure side — but only when honestly scored, or it
   manufactures phantom passes at a measured rate of 2 in 3.
   ⚠ **Weakened by v0.4.9's measurement.** More compute is *worth less* than this thesis implies:
   the `≥ 300 s` band spends half the sweep's CPU for a **1.0%** honest pass rate, and 93.1% of all
   passes already land under 30 s. Compute buys the tail very little. Re-read thesis 3 as "compute
   is *some* of the accuracy work, and the cheapest part of it is already bought."
4. **Canonicality is the larger axis and is generator-independent** — 520 molecules / 10.40 points vs
   340 failures / 6.80 points.
5. **Information-adding tokens go last.** P1/P2/P3 convert silent collapses into loud failures;
   promoting them lowers the headline until the generator reproduces them.

---

## The ladder

| ver | theme | headline prediction |
|---|---|---|
| **v0.4.7** | close-out: land the 5 finished swimlanes; tag `v0.4.7`; bump `pyproject` 0.4.6 → 0.4.7 | flat |
| **v0.4.8** | **the honest number** — promote `OIN_INDEP_SCORE` to *the* score; 5k honest sweep; + the atom-count gate | **DOWN ~5.7 pts, planned** |
| **v0.4.9** ✅ | **speed becomes measurable** — `OIN_ENFORCE_BUDGET` (default OFF); 328-molecule stratified benchmark frozen. **Refuted its own premise:** `elapsed_s` is a sum, the harness already enforces to ε ≈ 0.2 s | pass flat ✓; benchmark reproduces to **0.28%** |
| **v0.4.10** ✅ | **cost per attempt**, byte-identical. Deleted the discarded `.index()` scan (**`CAHQEJ` −32.9%**, `FOSNEI` +0.3% nil) and added `OIN_MEMO_CIP_REPARSE` (**`VAFMIA` −86.7%**, `CAHQEJ` −2.4%), both **byte-identical on ARM 1 and ARM 2**. ⚠ **Found its own arbiter broken:** ARM 1 had been non-runnable since `dd51a515` | pass FLAT ✓ by construction; **every speed number bimodal by molecule** |
| **v0.4.11** ✅ | **`slot_renumber`** (496 / 9.92 pts) — built the within-fragment donor fold v0.4.5 specified. 🔴 **REFUTED IT.** It works (+7.86 pts, 393 molecules, one direction) and **collapses enantiomers: 221 of its own 393 gains (56.2%), 19/250 on a uniform draw**, 18 oracle-confirmed chiral. Ships **default OFF**. Lane 1 classified all 496: `diff_occupancy` **0**, and 90 of the 118 `distinct_donors_LOCAL` are frozen-resonance-form artifacts ⇒ **v0.4.14, not here** | **FLAT** — the lever that would raise it is unsafe |
| **v0.4.12** | **honest acceptance** — harden `OIN_ATTACH_CHECK`; the `OIN_ETA_EARLY_EXIT` corpus A/B | biggest tail move; **pass must not move — that is the gate** |
| **v0.4.13** | **harness false negatives** — `PREFILTER_VETO` prevalence; the MEDZUR class | pass **UP** |
| **v0.4.14** | **`rdkit_canonical`** (114 / 2.28 pts) + winding residue — **RE-SIZED by v0.4.11 to ~4.08 pts**: 90 of the 118 residual `slot_renumber` molecules (**1.80 pts**) are the same ligand-BODY problem, a frozen resonance form (acac written ketone/enol) that `CanonicalRankAtoms` reads as two inequivalent donors | byte_exact up ~**4.1** pts |
| **v0.4.15** | **the 57 notation molecules** — build the per-case oracle | knowledge, not points |
| **v0.4.16** | **information-adding tokens** — P1 `\|mc:±\|`, P2 `\|ax:±\|`, P3 `[N@]` | **DOWN, then recovers** |
| **v0.4.17** | **`structural`** (417 / **8.34 pts**) — **RE-LABELLED**: generated structures that re-perceive with *different coordination than the OIN claims*. A generator-capability problem, not "wrong isomer"; it does not belong beside `facmer_divergent` | byte_exact **UP**, bounded by what the generator can assemble |
| **v0.4.18** | **generator capability floor** — boron assembly, `NON` geometry, 28 produced-nothing. **Encode floor R3** (15 / 0.30 pts) folded in here opportunistically | likely a documented limitation |

### LADDER DECISION 2026-07-27 — accepted, by the project owner

```
Change:   v0.4.11  was  encode floor R3 (15 mol / 0.30 pts)
                   is   slot_renumber   (496 mol / 9.92 pts)
          v0.4.14-v0.4.18 shift up one; `structural` pulled out of v0.4.16 and re-filed
          at v0.4.17, re-labelled as generator capability; encode floor R3 folded into
          v0.4.18 as opportunistic work.

Because:  the ladder pointed its NEXT release at the SMALLEST block in the gap while
          18.26 of 27.54 points sat four and six releases out. 0.30 points is a rounding
          error against a 27.54-point gap, and encode floor R3 depends on nothing and
          unblocks nothing, so no prerequisite argument protects its position. Its real
          merit -- ARM 1 of the byte-identity gate takes ~15 min, almost all of it
          fixtures hitting `_resonance_candidates_isolated` at 120 CPU-s per fork -- is a
          developer-convenience argument, not a roadmap-points one.

Not changed, and why:
          v0.4.12 (honest acceptance) and v0.4.13 (harness false negatives) stay early.
          They are MEASUREMENT PREREQUISITES: shipping a 9-point canonicality fix before
          the harness's false-negative rate is known risks attributing an artefact to the
          fix. That argument covers 12 and 13. It does not cover a 0.30-point release.

Also corrected: the ladder carried STALE counts -- `slot_renumber (459 / 9.18)` and
          `rdkit_canonical (61 / 1.22)` predate the v0.4.8 honest baseline. The honest
          figures are 496 / 9.92 and 114 / 2.28, reproduced 2026-07-27 from
          `results-v0.4.8-honest/bucket_report_honest.md`.
```


v0.4.8–v0.4.10 are written as full charter sets. v0.4.11–v0.4.18 are **SKETCH** and are promoted to
full one at a time by the close-out ritual below.

---

## Self-propagation

Every release ends by running `spec/handoffs/roadmap-100-100/CLOSEOUT.md`:

1. land the lanes · 2. re-measure and freeze the bucket report · 3. **diff predicted vs actual — a
miss is a deliverable** · 4. write the git-durable `docs/*.md` · 5. tag, bump, CHANGELOG ·
**6. generate `v0.4.<N+1>/` and promote it from SKETCH to full using the numbers just measured** ·
7. prune worktrees.

Step 6 runs `spec/handoffs/roadmap-100-100/PROMPT-next-release.md`, which requires the next session
to re-derive the gap decomposition, **re-check whether the ladder is still right**, and say
explicitly what the new numbers refuted. **v0.4.19 and beyond never need to be written by hand.**

---

## Standing traps

- **A sample that only exercises the common case confirms whatever you already believed.** On
  2026-07-26 a 10-molecule boron sample read "0/10 generate" and was used to propose a blanket
  fast-fail; the full 33 read **2/33** — `RAWJEG` (LIN, 2 slots) and `ULODUU` (TET, 4 slots) do
  assemble — and **the proposal is refuted**, because a blanket gate costs 2 real passes and no clean
  discriminator exists.
- **Circular measurement.** If the arm that scores an A/B uses the predicate the change accepts on,
  it cannot detect what the change costs. `OIN_ACCEPT_SCORED` read "18/22 both arms, zero
  regressions"; the honest arm read **15/20 → 7/20, 8 regressions, 0 fixes**.
- **Canonicalization can destroy what it canonicalizes.** Sorting symmetry-equivalent axes by sign
  made the axial token reflection-invariant and silently killed the chirality it encoded; every guard
  written against the single-axis fixture passed. **Prove the descriptor still flips for the mirror,
  on a multi-center case.**
- 🔴 **`byte_exact` CAN BE RAISED BY DELETING INFORMATION — and the comparison key will agree.**
  v0.4.11's donor fold moved **393 molecules into `byte_exact`, 0 in any other direction**, changed
  the key on **0 of 992** strings, and passed both gate arms — while **collapsing enantiomers in 221
  of those same 393 gains (56.2%)**. The metric and the key are *both* blind to reflection, because
  `compare._parse_vertex_colors` folds that axis deliberately. A one-directional transition matrix
  is therefore **not** evidence of safety for anything touching canonicalization.
  **Mirror-audit every canonicality lever (`tools/mirror_audit_donor_fold.py`) before quoting its
  points**, and confirm chirality independently with `tools/injectivity/oracle.py`. Two independent
  cohorts found it at comparable rates (uniform 19/250 = 7.6%, stratified 31/300 = 10.3%).
- **A partial run is not a result, and a progress line is not a tally.** During v0.4.11 the
  stratified audit was read mid-run as "0 regressions in the first 200" and written up as evidence
  that the cohort was the *wrong sample*. The tool prints a verdict only every 50th molecule, so
  four clean progress lines had been mistaken for 200 clean molecules; run to completion the same
  cohort reported **31 collapses**. The "wrong stratum" conclusion was an artifact of reading a
  sampled log as a census. **Wait for the summary line, or count with `--verbose`.**
- **A fragment's automorphism says nothing about the parity of the vertex permutation it induces.**
  Two donors interchangeable in the *isolated ligand graph* can sit at vertices whose exchange is an
  **improper** operation on the coordination sphere. This is why folding is restricted to proper
  rotations; narrowing a fold's scope *within a fragment* is not a substitute for a parity check.
- **`OIN_EMIT_AXIAL`'s promoting evidence is stale** — measured with `OIN_CANONICAL_PERCEPTION` OFF,
  now default-ON; YESKOZ's hindered-axis count goes 2 → 1 under it.
- **P3's obvious fix is measured-wrong** — stamping `[N@]` after the sanitize moves the canonical
  write order, and `@`/`@@` is a parity relative to that order; RIFGUJ's pseudo-asymmetric ring
  carbons then flip between a structure and its mirror.
- **vdW clash count is not a sufficient quality proxy** — `POVPIA_comp_0` goes 16 clashes → 0 while
  re-perception shows a detached hydrogen and a C–N read as C=N.
- **"Unset means off" is a trap.** `os.environ.get("OIN_EMIT_AXIAL")` returns `"0"`, a non-empty
  string, and *enables* the lever. Always go through `oin/levers.py::lever_enabled`.
- **The shared checkout moves under you** — one `.git/index` across ~33 worktrees.
- 🔴 **`metrics.elapsed_s` is a SUM, not a duration** (v0.4.9). Up to three separately SIGKILLed
  harness attempts are added into one field, so the ceiling is `3 × mol_timeout`. The "759.9 s
  against a 300 s budget" figure that chartered a whole release is `300 + 300 + 160`; split by
  `tier_passed`, all 4658 single-attempt rows finish within **0.2 s** of their cap. This is the
  **second** trap in the same field, and it survived one audit because everyone was busy
  remembering the first one (it is nested).
- 🔴 **A string-equality count is not a pass count** (v0.4.9). Comparing `smiles_1` to the honest
  round-trip string gives 3631 passes; the bucket report gives 3623, because it applies the
  `status` gate first. The eight are exactly the **atom-count gate** population — byte-identical
  strings with atoms missing. Use `tools/roundtrip_bucket_report.py`, never an ad-hoc
  `honest_class.endswith("->byte")`.
- 🔴 **A gate can silently run on the wrong interpreter** (v0.4.9). `gate_v047.sh` globbed
  `$(dirname $REPO)/*/.venv/bin/python` because worktrees have no `.venv`; from a worktree it
  selected an unrelated project's venv, then one with rdkit **2025.09.2** against the pinned
  **2025.9.3**. A byte-identity gate on a different rdkit reports MISMATCHes that read as code
  regressions. Resolve via `git rev-parse --git-common-dir` and refuse on version drift.
- 🔴 **`grep -c '[p]rocess_name'` matches your own shell** (v0.4.9). A wait-loop built on it never
  exits — twice in one session, once in a leftover watcher from an earlier session that had been
  spinning for hours. Match on something the watcher's own command line cannot contain.
- **Profile before bounding.** v0.4.9's charter named the unbounded CBC solve and the 48–57 s
  `accept_fn` re-encode as the cost sinks. Measured: **2.1%** and **0.8%**. The sink is
  `embed.get_embedding` (61.5 s of *self* time in an 82.4 s generation). Bounding either suspect
  would have measured as "no change" and shipped nothing while looking like a fix.

---

## Sources

`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md` · `docs/agentic-notes/v0.4.6/METRIC_FALSE_POSITIVES.md` · `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` §10 ·
`docs/agentic-notes/v0.4.5-retrospective/` · `src/oinsmiles/oin/levers.py` ·
and on `swimlane/v047-*` until v0.4.7 lands: `docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md`,
`docs/agentic-notes/v0.4.7/ATTACH_CHECK_v0.4.7.md`, `docs/agentic-notes/v0.4.7/ENCODE_FLOOR_v0.4.7.md`, `docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md`,
`docs/agentic-notes/v0.4.7/BORON_GEN_CEILING_v0.4.7.md`
