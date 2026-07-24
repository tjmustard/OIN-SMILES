# OIN-Direct Assembly — Validation & Decision (v0.4.4)

Record of the hybrid experiment that made **OIN-direct assembly the default** generation path
(replacing the lossy m-SMILES bridge) and rejected the winding-construction follow-on.
Raw per-arm data + drivers are in the gitignored `tmCAT-tmPHOTO_xyz_dataset/results-v0.4.4-promote-ab/`.

## Motivation

The generator built its internal `MetalComplex` by serialising the parsed OIN to a MetalloGen
`metal|lig|…|geo` m-SMILES string and re-parsing it — a redundant parse→serialize→parse round
trip through a format that is **lossy** (cannot express eta-ring winding; drops metal `@SPn`
chirality). Goal: build the complex directly from `ParsedOIN` so OIN-format changes reach 3D
generation without a translation layer, and winding + metal chirality survive.

## Experiments (standalone per-molecule driver, `seed=42`, fac/mer-key buckets)

**Step 1 — direct assembly + DG embed (`OIN_DIRECT_DG`).** Build the `MetalComplex` from
`ParsedOIN` (`om.get_om_from_parsed`), but run it through the SAME DG embed + winding-search +
early-exit body as the m-SMILES path. A/B on a 38-molecule stratified sample:

| arm | byte-exact | key-match |
|---|---|---|
| current (m-SMILES) | 60.5% | 73.7% |
| direct_dg | 60.5% | 73.7% |

→ **Neutral: 0 regressions, 0 gains, identical buckets.** Direct assembly is outcome-equivalent
to the m-SMILES path; it also *preserves* the two eta molecules SL2's rigid path regressed.

**Step 2 — deterministic winding construction (post-embed in-place ring flip).** With the
winding carried by direct assembly, flip each eta ring to its target winding after the DG embed
(reusing `_place_haptic`'s `signed_circulation` math). A/B on a 34-molecule eta sample (22
current failures + 12 passing-eta guards):

| arm | byte-exact | hard_fail | median |
|---|---|---|---|
| current | 47.1% | 1 | 34s |
| step2 (flip) | 44.1% | 7 | 37s |

→ **Clean NEGATIVE: 0 accuracy gains, +1 regression, +6 timeouts, slightly slower.** The eta
round-trip failures are **not** winding-search misses — the default's search + early-exit
already gets the ring face right, so constructing winding recovers nothing, and the per-conformer
flip cost tips slow eta molecules over the timeout. This is the **third** confirmation of
"selection beats construction" (after SL2 rigid rebuild and SL3 greedy placement, both negative).

**Confirmation — direct as the default.** Reverted Step 2; flipped `OIN_DIRECT_DG` to default-on
(m-SMILES retained as fallback). Re-A/B on the 38-molecule sample:

| arm | byte-exact | buckets (byte/key/facmer/struct/hard) |
|---|---|---|
| `OIN_DIRECT_DG=0` (old m-SMILES default) | 60.5% | 23 / 5 / 2 / 3 / 5 |
| default (new direct) | 60.5% | 23 / 5 / 2 / 3 / 5 |

→ **Identical: 0 regressions, 0 gains, 0 bucket changes.** Full unit suite **551 OK / 3 skip**
under the new default (goldens byte-identical via the direct path).

## Decision

- **PROMOTE OIN-direct assembly to default-ON** (`OIN_DIRECT_DG=0` / `ff_params["direct_dg"]=False`
  to opt out). The m-SMILES bridge (`convert_parsed_to_msmiles`) is retained **only as a
  fallback** — used if direct assembly raises on an edge case — so a rare case can't hard-fail.
  Rationale: removes the lossy translation layer and lets OIN-format changes affect 3D generation
  directly, at zero accuracy cost. Metal `@SPn` chirality now rides through the representation
  (inert until the encoder can reproduce it).
- **DROP Step 2 (winding construction).** Net-harmful; the residual eta failures are structural,
  not winding. Reverted entirely.
- **SL2 `oin_direct` (rigid haptic winding) stays opt-in** — regresses eta.

Guarded by `tests/unit/test_oin_direct_assembly.py::test_direct_dg_on_by_default`.
See `docs/GENERATION_PIPELINE.md` for the full default pipeline.
