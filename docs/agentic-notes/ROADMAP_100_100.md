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
| **A — accuracy** | `byte_exact`, scored by **independent** re-perception of the generated XYZ | 82.80% *(dishonest)* | **100%** |
| **B — speed** | `max(elapsed_s)` against an **enforced** budget | 759.9 s; **994/5000 = 19.88%** over 30 s | **< 30 s, p100** |

⚠ `metrics.elapsed_s` is **nested**; read from the top level it silently yields `0`.

## The gap — `100 − 82.80 = 17.20` points

| block | n | pts | nature | release |
|---|---:|---:|---|---|
| `key_equal` → `slot_renumber` | 459 | 9.18 | canonicality, encoder-side | v0.4.14 |
| `key_equal` → `rdkit_canonical` | 61 | 1.22 | canonicality, encoder-side | v0.4.15 |
| `hard_fail` | 315 | 6.30 | compute (mostly) | v0.4.9 – v0.4.13 |
| `encode_fail` | 15 | 0.30 | encoder coverage | v0.4.11 |
| `structural` | 9 | 0.18 | wrong isomer | v0.4.16 |
| `facmer_divergent` | 1 | 0.02 | wrong isomer | v0.4.16 |
| **sum** | **860** | **17.20** ✓ | | |

Plus, outside the table: **the ~5.7-point metric over-statement**, owned by v0.4.9's predecessor.

Runtime, recomputed from the primary reports: overall n = 5000, median **7.19 s**, `> 30 s` **994
(19.88%)**, max **759.9 s**; **eta** n = 1146, median 24.08 s, `> 30 s` 528 (46.1%) — **53.1% of the
whole tail**; non-eta n = 3854, median 6.18 s, `> 30 s` 466 (12.1%).

---

## Ordering thesis

1. **Instrument before engine.** The headline goes **down once, on purpose**, inside its own release
   boundary — otherwise it is indistinguishable from a regression, a confound this project has been
   caught by twice.
2. **Bound before optimize.** A p100 target against an unbounded budget is not a target.
3. **Compute is the accuracy work** for the failure side — but only when honestly scored, or it
   manufactures phantom passes at a measured rate of 2 in 3.
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
| **v0.4.9** | **speed becomes measurable** — enforce the budget; freeze a stratified runtime benchmark | pass flat; tail bounded |
| **v0.4.10** | **cost per attempt**, byte-identical — SVD in `_finalize_positions`; redundant per-attempt work | pass flat *by construction*; median down |
| **v0.4.11** | **encode floor R3** — memoize the forked-resonance timeout verdict per ligand graph | `encode_fail` down |
| **v0.4.12** | **honest acceptance** — harden `OIN_ATTACH_CHECK`; the `OIN_ETA_EARLY_EXIT` corpus A/B | biggest tail move; **pass must not move — that is the gate** |
| **v0.4.13** | **harness false negatives** — `PREFILTER_VETO` prevalence; the MEDZUR class | pass **UP** |
| **v0.4.14** | **`slot_renumber`** (459 / 9.18 pts) | byte_exact **UP ~9 pts** |
| **v0.4.15** | **`rdkit_canonical`** (61 / 1.22 pts) + winding residue | byte_exact up ~1.2 pts |
| **v0.4.16** | **the 57 notation molecules** — build the per-case oracle | knowledge, not points |
| **v0.4.17** | **information-adding tokens** — P1 `\|mc:±\|`, P2 `\|ax:±\|`, P3 `[N@]` | **DOWN, then recovers** |
| **v0.4.18** | **generator capability floor** — boron assembly, `NON` geometry, 28 produced-nothing | likely a documented limitation |

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

---

## Sources

`docs/agentic-notes/v0.4.6/SWEEP_v0.4.6_5K.md` · `docs/agentic-notes/v0.4.6/METRIC_FALSE_POSITIVES.md` · `docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md` §10 ·
`docs/agentic-notes/v0.4.5-retrospective/` · `src/oinsmiles/oin/levers.py` ·
and on `swimlane/v047-*` until v0.4.7 lands: `docs/agentic-notes/v0.4.7/ACCEPT_SCORED_v0.4.7.md`,
`docs/agentic-notes/v0.4.7/ATTACH_CHECK_v0.4.7.md`, `docs/agentic-notes/v0.4.7/ENCODE_FLOOR_v0.4.7.md`, `docs/agentic-notes/v0.4.7/COHORT_v0.4.7.md`,
`docs/agentic-notes/v0.4.7/BORON_GEN_CEILING_v0.4.7.md`
