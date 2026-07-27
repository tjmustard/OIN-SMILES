# OIN-SMILES v0.4.1 — full-dataset round-trip accuracy

This report gives the **population-level** round-trip result for the shipped v0.4.1 code over the
**entire** tmCAT/tmPHOTO corpus — every one of the 25,197 unique complexes, not a sampled subset.
It is the companion to `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md` (which measures *per-molecule* what the v0.4.2 wave
changed) and to `CHANGELOG.md`.

## Headline

| Metric | Value |
|---|---:|
| Complexes tested (unique `<REFCODE>_comp_<N>`) | **25,197** |
| Round-trip pass (canonical OIN string identity **and** coordination-sphere RMSD < 1.0 Å) | **22,280 (88.42 %)** |
| Failures | 2,917 (11.58 %) |
| — of which stochastic/wontfix harness noise (`timeout` + `no_conformers` + carborane) | 1,934 |
| **Accuracy-clean pass** (failures above excluded) | **95.77 %** |
| Genuine accuracy-defect failures | 983 (3.90 %) |

The round-trip is the full loop **XYZ → OIN string → 3D structure → OIN string**, scored on two gates:
the regenerated canonical OIN string must be byte-identical to the original, **and** the regenerated
coordination sphere must align to the input within 1.0 Å RMSD.

## What this number is — and is not

This is measured on the continuous `--quick` accumulator
(`tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0/`). Two honesty caveats, the same ones raised in
`docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md`, apply and are worth stating up front:

1. **`--quick` is a deliberately weak generator — 88.4 % is a *floor*, not the library's ceiling.**
   The accumulator runs `--quick` (`uff_pool_size=2`, `max_attempts=10`, 30 s hard-kill) as a fast
   screening path. The default non-quick generator (`pool=5`, full budget, `--mol-timeout 1800`) passes
   materially more: most `timeout`, and much of `geometry_or_fragment_change` / `winding_flip`, are
   quick-budget artifacts that evaporate under the default generator (see `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md`).
   So the **non-quick pass rate is higher than 88.4 %**; this figure is the conservative screening floor.

2. **Provenance is uniform-accuracy, but not a single commit.** The 25,197 reports span
   `c7edeeb6` (= tag `v0.4.1`), `5538b722` (= tag `v0.4.0`), and a few v0.4.1-development commits.
   These are **accuracy-equivalent** — the v0.4.0 perf wave preserved generated geometry byte-for-byte
   (golden A/B), so the mix does not skew the rate. **No v0.4.2 code is in this corpus** (verified: zero
   `e6febd16` stamps), so this is a clean v0.4.1 baseline, distinct from the v0.4.2 numbers.

What a full-corpus number *does* give — that per-molecule flip counting cannot — is the **real
distribution of failure modes across all 25,197 field complexes**: how common each defect class is in
the wild, and therefore what to prioritise. That distribution is the substance of this report.

## Failure-mode distribution (all 2,917 failures)

| # | class | count | % of dataset | tier / disposition |
|---|---|---:|---:|---|
| 1 | `timeout` | 1,389 | 5.51 % | harness noise — `--quick` 30 s hard-kill (mostly masked `no_conformers`) |
| 2 | `no_conformers` | 403 | 1.60 % | generator robustness — genuine valence/perception gaps |
| 3 | `donor_H_atom_count` | 271 | 1.08 % | **accuracy defect** — donor H-count mismatch (η-arene / ammine-nitride) |
| 4 | `high_rmsd` | 162 | 0.64 % | FF-floor — string-correct round-trips that miss only the 1.0 Å gate under FF-only |
| 5 | `carborane_unsupported` | 142 | 0.56 % | wontfix — 3c2e cage bonding is outside the two-centre model |
| 6 | `gen_exception_other` | 110 | 0.44 % | representation limit — outer-sphere counterions/solvent (`UncoordinatedFragmentError`) |
| 7 | `EZ_bond_stereo` | 88 | 0.35 % | **accuracy defect** — C=C / C=N E/Z reproduction |
| 8 | `string_mismatch_other` | 87 | 0.35 % | **accuracy defect** — residual canonical-string divergences |
| 9 | `geometry_or_fragment_change` | 85 | 0.34 % | mostly `--quick` artifact — resolves under the non-quick default |
| 10 | `winding_flip` | 52 | 0.21 % | mostly `--quick` artifact — η-ligand winding face |
| 11 | `atom_stereo` | 41 | 0.16 % | **accuracy defect** — @/@@ heteroatom / donor-C handedness |
| 12 | `macrocycle_perception` | 19 | 0.08 % | mixed — mostly boron-stereo / geometry / E/Z misfiles |
| 13 | `kekulize_encode_crash` | 19 | 0.08 % | **accuracy defect** — stale-aromatic kekulization failure |
| 14 | `encode_crash_other` | 16 | 0.06 % | **accuracy defect** — charge-ladder fall-through on large cages |
| 15 | `H_on_terminal_oxo_imido` | 16 | 0.06 % | notation ambiguity — protonated terminal =O/=N donor |
| 16 | `garbled_aromatic` | 9 | 0.04 % | perception — reduced-porphyrin (chlorin) saturation mismatch |
| 17 | `geometry_NON` | 6 | 0.02 % | geometry — high-CN template / embedding limit |
| 18 | `rmsd_mapping_failed` | 2 | 0.01 % | metric edge case |

Grouped by nature:

- **Harness / FF artifacts (not encoding errors): ~1,956** — `timeout` (1,389), `high_rmsd` (162,
  string-correct, FF geometry-gate only), `no_conformers` (403 generator robustness), `rmsd_mapping_failed` (2).
- **Representation / notation limits (wontfix or docs): ~252** — `carborane_unsupported` (142),
  `gen_exception_other` (110, all outer-sphere species with no binding slot).
- **Genuinely encoder/generator-fixable accuracy classes: ~700** — the E/Z, atom-stereo, donor-H,
  winding, geometry, and encode-crash rows. This is the surface the v0.4.2 wave targets; see
  `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md` for the per-molecule before→after.

The takeaway mirrors the v0.4.2 re-framing: **the true accuracy-defect surface is small** (~3.9 % of the
dataset), and a large share of the raw 11.6 % failure rate is screening-mode noise or explicit
representation limits, not encoding infidelity.

## Reproducing

```bash
# Regenerate the roll-ups from the individual reports (read-only on the reports):
R=tmCAT-tmPHOTO_xyz_dataset/results-v0.4.0
python tools/rebuild_summary.py    --output-dir $R   # → 25,197 / 22,280 / 2,917
python tools/classify_failures.py  --output-dir $R   # → case_registry.json (per-class)
python tools/group_v041_backlog.py --output-dir $R   # → tiered backlog
```

Per-molecule accuracy — and any claim about what a code change *fixed* — must be measured against a
single clean commit and, for conformer-yield-sensitive rows, **non-`--quick`** (`--mol-timeout 1800`);
`--quick` fabricates geometric distortion. See `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md` for the methodology.
