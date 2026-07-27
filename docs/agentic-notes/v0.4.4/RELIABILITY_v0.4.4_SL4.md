# SL4 reliability — acceptance report (v0.4.4)

Swimlane `swimlane/v044-reliability`. Target: shrink the v0.4.4 BASELINE `hard_fail` cohort
(**306 molecules / 4.55%**; the honest post-SL0 baseline is 93.51%) without regressing the passing
middle. Decisions: **demote RMSD to a diagnostic**, **gated no-progress cutoff (OFF by default)**,
**FF-only** measurement.

## Levers shipped

1. **RMSD demoted to a diagnostic** (`tools/test_dataset_roundtrip.py::_attempt_generation`). The
   canonical-key string match is the lossless-hash pass contract; coordination-sphere RMSD is only
   ~0.22-correlated with geometric quality (`FALSIFICATION_v0.4.3_ELIMINATION` §1.4). A string-exact
   round-trip that only exceeds `RMSD_GATE` is now a **success** carrying an `rmsd_over_gate`
   diagnostic, not a failure. `rmsd is None` (unmappable sphere) still fails — a metric error, not a
   demonstrated geometry pass.
2. **`clash_count` diagnostic** (via `vdw_clash_count` on the generated coords). The "no-clash" half
   of the contract is enforced upstream by the generator's vdW acceptance gate; a hard clash gate in
   the harness is infeasible while `best_rejected_returned` ships the least-bad conformer ~86% of the
   time, so clash is **measured, not gated**.
3. **Gated no-acceptance-progress cutoff** (`generator3d/__init__.py`, opt-in
   `ff_params["embed_no_progress_attempts"]` / `OIN_EMBED_NO_PROGRESS`). Stops the attempt loop after
   N consecutive gate-rejected embeds (the OSIHUU pattern). **OFF by default → success path
   byte-identical.** This is a latency lever for a future A/B; it contributes **0** to the accuracy
   number here by design.

## Method

FF-only (xtb forced unavailable to keep the comparison honest against the frozen FF baseline and to
avoid the latent multiplicity defect a live xtb activates), `seed=42`, `--mol-timeout 1700`, 6
shards, 306 molecules. Raw reports + reproduction scripts:
`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.4-sl4/` (gitignored). Classified with
`tools/triage_hard_fails.py` on the SL0 key.

## Result — cohort 306 → 255 (−51, −16.7%)

| outcome | count | eta | non-eta |
|---|---:|---:|---:|
| **reclaimed (now key-exact success)** | **51** | 37 | 14 |
| still failed | 255 | 70 | 185 |

**Reclaim attribution (honest):**

| reason | count | what it is |
|---|---:|---|
| `rmsd_demote` | 28 | string-exact structures the old harness rejected on RMSD alone — **causally this lever** |
| `passed_outright` | 23 | stored false-failures that pass the RMSD gate cleanly on a fresh FF run — re-run drift (determinism / rdkit-version between the capstone commit and now), **not** this lever, but they show the frozen baseline over-counted failures |

Reclaimed structures are clean: **46/51 (90.2%) clash-free**, mean `clash_count` 0.37.

**Still-failed populations (disjoint) — outside SL4's levers:**

| population | count | eta | non-eta |
|---|---:|---:|---:|
| `atom_count` | 113 | 24 | 89 |
| `no_conformer` | 109 | 35 | 74 |
| `gen_exception` | 25 | 6 | 19 |
| `string_mismatch` | 5 | 5 | 0 |
| `timeout` | 3 | 0 | 3 |

`atom_count` (wrong H count) and `no_conformer` (unassemblable / budget-limited pool) dominate — these
are generation-quality failures (SL2/SL3 + the om.py/adapter path), not reliability-gate artifacts.
`string_mismatch` (5, all eta) are genuine positional/winding round-trip errors (SL2/SL3 territory).

## Eta overlap (SL2)

Reclaims are eta-heavy (37/51), and **70 of the 255 still-failed are eta** — the timeout/winding
cohort SL2's oin-direct winding targets. Per the handoff, re-measure this eta residual after SL2
lands; do not attribute it to SL4.

## Acceptance gate

1. **Cohort shrinks — yes**, 306 → 255, with attribution above (28 lever / 23 drift).
2. **No regression on the passing middle — yes, by construction.** The no-progress cutoff is OFF by
   default (success path byte-identical, pinned by `test_reliability_budget.py`); the RMSD demote only
   converts failures to passes, never the reverse. The 10-molecule pre-sweep confirmation showed 4/4
   passing molecules unchanged.
3. **Determinism preserved** (`seed=42`; unit test).
4. **Full unit suite green — 504 OK / 3 skip** (496 floor + 8 new); `ruff` clean.

## Follow-ups (out of this swimlane)

- The `rmsd is None` (unmappable-sphere) rows still fail; a future pass could demote these too once
  the mapping-failure cohort is characterized.
- A/B the no-progress cutoff (ON vs OFF) on the `no_conformer` cohort to measure latency savings
  without regressing the passing middle, then decide on promotion.
- Raising the internal embed budget above 300 s to test the >300 s healthy-but-slow hypothesis
  (kept at 300 s here to match the frozen baseline).
