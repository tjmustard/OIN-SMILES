# v0.4.4 — Promote A/B (integration decision)

**Run 2026-07-23**, `main` `30bb5471` (all 6 swimlanes landed, all accuracy levers gated OFF).
Measurement-only; this folder is gitignored (durable copy of the `/tmp` scratchpad run).

## Question

v0.4.4 landed three generation-accuracy levers gated OFF: SL1 `early_exit`
(`OIN_EARLY_EXIT`), SL2 `oin_direct` (`OIN_DIRECT_ASSEMBLY`), SL3 `greedy`
(`OIN_GREEDY_PLACEMENT`) (+ SL1 `stretched_bond`). Which, if any, ship default-ON?

## Method

Standalone per-molecule driver (`ab_driver.py`), one deterministic conformer (`seed=42`,
`optimizer="ff"`) per molecule per arm; re-encode via `get_oin_string(contract_mol)` (the
BASELINE.md / harness metric); bucket with SL0's `roundtrip_bucket_report.classify`. Each arm
= the same driver with a different env (levers are env-controllable); vdW gate stays ON (shipped
default) in all arms. NOT the two-pass harness (its g-xTB recovery would contaminate the FF
arms — A5's lesson).

**Sample: 38 molecules, stratified by fac/mer bucket** (`sample.txt`) — deliberately
worst-cohort-biased (44.7% baseline byte-exact vs 81% population) so the levers have room to
act: 6 byte-exact non-eta + 6 byte-exact eta (regression guards) + 4 key_equal + **all 12
facmer_divergent** + 6 structural + 4 hard_fail.

**Promote rule (fixed before results):** ship default-ON only if a lever **raises byte-exact /
key-match with ZERO regressions** (no baseline byte-exact molecule breaks), or is
quality-neutral with a clear **speed** win. Anything that regresses stays opt-in.

## Results

| arm | byte_exact | key_equal | facmer_div | structural | hard_fail | byte% | key% | med s | sum s |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| A_baseline    | 17 | 4 | 7 | 3 | 7 | 44.7 | 55.3 | 53 | 3739 |
| **B_early_exit** | **23** | 5 | 2 | 3 | 5 | **60.5** | **73.7** | **10** | 2606 |
| C_oin_direct  | 15 | 7 | 6 | 3 | 7 | 39.5 | 57.9 | 29 | 3096 |
| D_greedy      | 18 | 6 | 4 | 3 | 7 | 47.4 | 63.2 | 51 | 3846 |
| E_all_on      | 16 | 8 | 3 | 3 | 8 | 42.1 | 63.2 | 23 | 2870 |

**Regressions / gains vs baseline** (regression = byte-exact → not):

| arm | byte regressions | byte gains | key-match gains |
|---|---|---|---|
| B_early_exit | **0** | 6 (FOKNOL, IVAWIL, KELQOI, MANKIX, SOXTEH, VOZKAY — all facmer) | 7 |
| C_oin_direct | 2 (ABEXOU, AKOXII — eta) | 0 | 3 |
| D_greedy | 1 (ABIROU — eta) | 2 (FOKNOL, MANKIX) | 3 |
| E_all_on | 3 (union of C+D) | 2 | 5 |

## Decision

| Lever | Decision | Evidence |
|---|---|---|
| **SL1 `early_exit`** | **PROMOTE to default-ON** (`OIN_EARLY_EXIT=0` to opt out) | byte-exact 44.7→60.5%, key-match 55.3→73.7%, **0 regressions** (incl. all 12 guard molecules), **~5× faster** (median 53→10s). Picks the conformer reproducing the requested fac/mer isomer; non-regressive by construction (falls through if no pool conformer matches). |
| **SL2 `oin_direct`** | **keep opt-in** | Net byte-exact DOWN (−2); **regressed 2 eta molecules** already byte-exact. Direct assembly perturbs winding on structures that were fine. Needs the eta regressions fixed first. |
| **SL3 `greedy`** | **keep opt-in** | +1 net but carries 1 eta regression and no speed benefit; its gains (FOKNOL, MANKIX) are already captured by early_exit. |
| **SL1 `stretched_bond`** | **keep opt-in** | Only exercised inside the (regressing) all-on arm; no clean isolated signal to justify default-on. |

**Takeaway (confirms A5's "selection beats placement"):** the *selection* lever wins cleanly on
both quality and speed; the *placement* levers (oin_direct, greedy) perturb already-correct
geometries and regress eta round-trips. Only `early_exit` is promoted.

## Reproduce

Driver + sample + per-arm reports in this folder. Re-run:
`bash run_arms.sh` (writes `<arm>/individual_reports/`), then `python analyze.py`.
