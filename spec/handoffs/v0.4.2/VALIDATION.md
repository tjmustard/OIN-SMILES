# v0.4.2 Capstone Validation — RESULTS (corrected)

**Date:** 2026-07-16 → 2026-07-17
**Branch under test:** `release/v0.4.2` @ `58bba7ad` · **Baseline:** `c7edeeb6` (= v0.4.1 = main)
**rdkit:** 2025.09.3 · **optimizer:** ff · **seed:** 42 (default)
**Verdict:** **GO for merge — zero production (quick-mode) regressions.** Merged on user go-ahead.

> This file was revised twice. The final verdict rests on the **corrected methodology** in §1; earlier
> drafts stated "GO via 90/90 identical" from a comparison that was silently v0.4.2-vs-v0.4.2 (see §1).
> The v0.4.0-style runbook is preserved below the `═══` separator.

## 1. Methodology correction (load-bearing)

**The generator that runs is set by the HARNESS's repo, NOT the venv.** `tools/test_dataset_roundtrip.py`
does `sys.path.append("../src")`, and the generation subprocess loads `oinsmiles` from the harness's
own `src/` — proven by a `GENERATE_FROM __file__` trace: the **main-repo** harness runs
`OIN-SMILES/src` (c7edeeb6) and the **worktree** harness runs `OIN-SMILES-v0.4.2/src` (v0.4.2).
(The parent-process `oinsmiles loaded from` debug line shows the editable venv path and is MISLEADING
for provenance; and the report `commit_id` stamp reflects the harness location, not the generator.)

Consequence: swapping the **venv** does NOT swap the generator. The correct c7edeeb6-vs-v0.4.2 A/B is
**main-repo harness vs worktree harness**. The capstone SWEEP (worktree harness) genuinely measured
v0.4.2 and is valid; the intermediate triage/fixconfirm/quickab (which swapped only the venv under the
worktree harness) compared v0.4.2 to itself and were discarded.

## 2. The sweep

6719-molecule non-quick sweep (worktree harness = v0.4.2), `--mol-timeout 1800`, coverage 6719/6719,
uniform provenance (`quick=False, ff, rdkit 2025.09.3`). Accumulator paused for the window (frozen).
5943/5964 baseline-quick passers still pass under v0.4.2 non-quick → 21 regressor candidates + 352
apparent flips into pass. These were then adjudicated with the CORRECT generator A/B (§3).

## 3. Corrected gate — main-repo(c7edeeb6) vs worktree(v0.4.2), 90 gate molecules

Ran the 21 candidate regressors + 69 code-class candidate fixes on the **c7edeeb6 generator**
(main-repo harness) non-quick, and cross-checked the 3 surviving regressors in a full 2×2
(generator × {quick, non-quick}):

| molecule | c7 quick | v0.4.2 quick | c7 non-quick | v0.4.2 non-quick | verdict |
|---|---|---|---|---|---|
| OWODEK | FAIL | FAIL | PASS (0.19) | FAIL | non-quick-only |
| TIPYEX | FAIL | FAIL | FAIL (1.14) | FAIL (1.14) | not a regression |
| YEPXID | FAIL | FAIL | PASS (0.72) | FAIL | non-quick-only |

- **Production (quick) regressions: 0.** All three fail on BOTH generators in quick — the mode the
  accumulator/production actually uses. v0.4.2 breaks nothing that c7edeeb6 generates in quick.
- **Non-quick-only regressions: 2** (OWODEK, YEPXID): c7edeeb6 non-quick lands a clean conformer that
  v0.4.2 non-quick misses; both still fail in quick, so not a production regression.
- **TIPYEX is not a regression:** run alone, c7edeeb6 non-quick also fails (1.1368, identical to
  v0.4.2). The earlier "c7 passes" was the **sharded** run selecting a luckier conformer (0.7931).
- **Proven v0.4.2 code fixes (c7 non-quick FAIL → v0.4.2 non-quick PASS): ~25** — the wave's real value.

## 4. Meta-finding: the metric is conformer-flaky at the margin

Round-trip pass/fail is sensitive to which conformer the embed selects, and that selection varies with
run context: TIPYEX flipped PASS↔FAIL between sharded and alone runs; the accumulator recorded
OWODEK/YEPXID as quick-PASS but isolated quick re-runs FAIL them. So BOTH the baseline (a single
`--quick` sample) and any single-run gate sit on a noisy substrate; marginal fixes/regressions are
largely conformer noise unless confirmed across multiple runs. A future hardening: pin conformer
selection on the eval path so the accuracy metric is reproducible, then re-baseline.

## 5. Verdict

**GO for merge.** Zero production (quick-mode) regressions; ~25 real non-quick fixes; the 2
non-quick-only regressions do not affect the production metric. Merged to `main` as a squash on the
user's explicit go-ahead, **unpushed**, no version bump. `spec/handoffs/v0.4.2` stripped from the main
commit (stays gitignored on main). Recommended follow-up: reproducible-conformer hardening (§4) before
trusting future single-run accuracy deltas.

═══════════════════════════════════════════════════════════════════════════════════════════════════
# ▶ START HERE — capstone validation (v0.4.2 round-trip accuracy wave)  [ORIGINAL RUNBOOK]

**Launch a fresh Claude Code session in the main checkout and hand it this file.** Run **only after
every fix phase (S1, S3, S5, S6a, S6b, S7, docs) has squash-merged into `release/v0.4.2`.** The
capstone **edits no source** — a capstone that edits the thing it measures is not a measurement. If
it finds a regression, it reports + opens a follow-up phase.

### 1 · Work in the staging worktree
```bash
cd /home/tjmustard/Documents/GitHub/OIN-SMILES-v0.4.2   # on release/v0.4.2, all phases landed
git log --oneline release/v0.4.2 ^c7edeeb6              # confirm every phase squash is present
uv run python -m unittest discover tests/unit           # MUST be green BEFORE measuring
```

### 3 · The gate — per-molecule set-inclusion, NOT a percentage
`{passes on release/v0.4.2} ⊇ {passes on c7edeeb6}`. Any c7edeeb6-passer that now fails is a named
blocker with root cause. Claim named flips per-class AND per-CN AND per-metal. **CRITICAL (learned
the hard way): the generator is set by the harness's repo `../src`, so the c7edeeb6 arm MUST use the
main-repo harness — swapping the venv is a no-op.** Confirm regressions in the mode production uses
(`--quick`), not only non-quick.

### 5 · Merge to main (only on GO, and on the user's go-ahead)
`git -C .../OIN-SMILES merge --squash release/v0.4.2`; strip `spec/handoffs/v0.4.2`; commit; leave
unpushed; version bump + tag only on a further go-ahead.
