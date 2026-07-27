# MetalloGen structure-distortion research (seed for v0.4.3)

Research report characterizing how physically distorted MetalloGen's generated 3D structures
are, to scope the **v0.4.3** structure-quality wave. It is measurement only — no MetalloGen
code was changed. It is reproducible via `tools/structure_distortion_report.py`; the
per-molecule metrics live in `<results>/distortion_metrics.json` and the auto-generated tables
in `<results>/distortion_report.md`.

Population: the **non-quick capstone sweep** `tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042/`
— **6,404 generated structures** (6,295 successes + 109 failed-with-geometry), each scored
against its matched real crystal-structure input. The `--quick` sweeps are excluded (they
fabricate distortion; see `docs/agentic-notes/v0.4.2/ACCURACY_v0.4.2.md`).

## TL;DR

- **Distortion is real, systemic, and steric — not bond-length.** **53%** of generated
  structures carry ≥1 van-der-Waals steric clash (mean ~2.4, worst 28) vs **5%** of the real
  inputs. Bond *lengths* are near-ideal; the problem is atoms and ligand bodies driven into each
  other, plus elevated bond-angle strain (6.3° vs 3.4°).
- **It is population-wide, not a failure-only problem.** Structures that pass the topological
  round-trip gate clash as much as the failures — the string/RMSD gates are blind to geometry,
  and the stored coordination-sphere RMSD (~0.33 Å) looks fine because the clashes are in the
  ligand bodies it never inspects.
- **MetalloGen's own clash gate is too permissive.** Its acceptance test forbids only atomic
  *fusion* (covalent-radius ratio < 0.6); **100%** of shipped structures pass it while half
  carry vdW clashes. This is the most direct v0.4.3 lever.
- **Worst strata (priority):** early, oxophilic, high-coordination metals — **Y (MPO 0.372,
  92% clashing), Sc (0.460, 97%), Zr, Ti, Ta, Nb, Hf** — and high-CN / low-symmetry geometries
  **PBP, SQA, TBP, SPY**. Cleanest: **Au (0.863), Hg, Ni, Zn** and **LIN (0.875)**.

## The MPO quality score

Each structure is scored 0–1 by an MPO (multi-parameter optimization) over three angles that
were **validated to discriminate** distortion (generated vs. real input). Each raw metric maps
to a Derringer–Suich desirability (minimize/maximize, bounds calibrated to the input
population); angle sub-scores are weighted means; the MPO is their weighted **geometric** mean
(a fully-failed angle strongly penalizes). Weights/bounds are `--config`-overridable.

| Score | median | p25 | p10 | components (weight) |
|---|---:|---:|---:|---|
| **MPO (overall)** | **0.698** | 0.595 | 0.462 | geometry 0.45 · graph 0.30 · reference 0.25 |
| geometry sub-score | 0.772 | 0.556 | 0.255 | vdW clashes, severe clashes, worst overlap, angle strain |
| graph sub-score | 0.556 | 0.500 | 0.389 | perceived-CN divergence, bond-count divergence, perception health |
| reference sub-score | 0.813 | 0.726 | 0.637 | coordination-sphere RMSD + full heavy-atom overlay divergence |

## Headline metrics — generated vs. real input

| Metric | Gen median | Gen p95 | Input median | Input p95 |
|---|---:|---:|---:|---:|
| vdW clashes (< 0.75·ΣR_vdw) | 1.00 | 10.00 | 0.00 | 1.00 |
| severe clashes (< 0.60) | 0.00 | 3.00 | 0.00 | 0.00 |
| worst non-bonded overlap ratio | 0.74 | 0.82 | 0.80 | 0.85 |
| bond-angle strain (deg) | 6.32 | 13.56 | 3.42 | 10.72 |

**53% of generated structures have ≥1 vdW steric clash, vs 5% of inputs.**

### Two metrics were computed but excluded from the MPO (they do not discriminate)

| Metric | Gen median | Input median | Verdict |
|---|---:|---:|---|
| bond-length deviation (frac) | 0.038 | 0.087 | **inverted** — MetalloGen places near-ideal bond lengths, so it scores *better* than real crystals. Bonds are not the distortion. |
| ligands-only UFF displacement (Å) | 1.01 | 0.96 | **confounded** — an all-single-bond, metal-free FF relaxes real inputs about as far (energy/atom is inverted). Measures bond-order-model mismatch, not distortion. N/A rate 11%. A valid "would it fail optimization" metric needs a metal-capable engine (xtb/MACE), absent in this env. |

## MPO by coordination geometry

| geometry | n | median MPO | % clashing | mean clashes |
|---|---:|---:|---:|---:|
| PBP | 37 | 0.515 | 89% | 6.24 |
| SQA | 17 | 0.575 | 71% | 5.88 |
| TBP | 280 | 0.640 | 68% | 3.75 |
| SPY | 372 | 0.661 | 67% | 3.31 |
| TPL | 311 | 0.664 | 61% | 2.76 |
| OCT | 1542 | 0.673 | 56% | 2.87 |
| TPY | 157 | 0.681 | 58% | 2.43 |
| TET | 1236 | 0.682 | 51% | 1.95 |
| SPL | 2070 | 0.726 | 50% | 1.96 |
| LIN | 382 | 0.875 | 27% | 0.63 |

## MPO by metal (worst 12 of 27)

| metal | n | median MPO | % clashing | mean clashes |
|---|---:|---:|---:|---:|
| Y | 71 | 0.372 | 92% | 9.14 |
| Sc | 30 | 0.460 | 97% | 7.63 |
| Ta | 19 | 0.502 | 68% | 6.74 |
| Nb | 17 | 0.561 | 82% | 5.35 |
| Hf | 50 | 0.586 | 68% | 3.68 |
| Ru | 686 | 0.624 | 73% | 3.97 |
| Cd | 30 | 0.625 | 53% | 1.80 |
| Os | 44 | 0.630 | 61% | 2.86 |
| Rh | 327 | 0.634 | 76% | 3.32 |
| Ir | 479 | 0.648 | 57% | 2.55 |
| Zr | 247 | 0.653 | 56% | 2.80 |
| Mo | 192 | 0.660 | 44% | 2.49 |

Cleanest metals: **Au 0.863, Hg 0.787, Ni 0.780, Mn 0.760, Zn 0.765, Co 0.739, Cr 0.736**. The
gradient tracks ionic radius / oxophilicity / coordination number — big early metals hold more,
bulkier ligands that the current placement crowds.

## Worst 25 by MPO — the v0.4.3 priority queue

| molecule | metal | geom | MPO | clashes | angle° | full-div Å | coord-RMSD |
|---|---|---|---:|---:|---:|---:|---:|
| NEXTIT_comp_0 | Zr | OCT | 0.085 | 9 | 20.1 | 1.13 | 0.53 |
| MEXXES_comp_0 | Zr | PBP | 0.085 | 11 | 21.8 | 0.90 | 0.62 |
| EFEHON_comp_0 | Y | OCT | 0.089 | 14 | 19.2 | 1.18 | 0.79 |
| OVIPUF_comp_0 | Sc | OCT | 0.153 | 15 | 17.6 | 1.00 | 0.69 |
| IVEMOL_comp_0 | Ru | TET | 0.160 | 6 | 17.4 | 1.03 | 0.57 |
| MELFOB_comp_0 | Zr | PBP | 0.161 | 10 | 17.5 | 0.99 | 0.42 |
| DIDSAM_comp_0 | Y | OCT | 0.162 | 12 | 17.1 | 1.13 | 0.78 |
| BUNTEK_comp_0 | Ti | OCT | 0.191 | 12 | 16.3 | 0.99 | 0.45 |
| QEXHEH_comp_0 | Zr | OCT | 0.192 | 6 | 16.0 | 1.08 | 0.65 |
| KOCCEN_comp_0 | Y | TBP | 0.196 | 12 | 15.9 | 1.12 | 0.57 |
| RILHOL_comp_0 | Y | OCT | 0.219 | 20 | 13.8 | 1.52 | 1.36 |
| PIYNIV_comp_0 | Y | OCT | 0.226 | 28 | 14.5 | 1.04 | 0.54 |

(Full ranked list of all 6,404 in `distortion_metrics.json`, sorted by `mpo`.) The tail is
dominated by Y/Zr/Ti/Sc octahedra and bipyramids with 6–28 clashes and 13–22° angle strain —
severe crowding, yet coordination-sphere RMSD stays ≤ ~0.8 Å, confirming the sphere metric
misses it.

## Root-cause signals & hypotheses for v0.4.3 (evidence-based, not yet fixes)

1. **The clash gate is too permissive.** `clean_geometry.py` accepts on `ratio_criteria=0.6`
   (covalent radii) + `atom_d_criteria=0.5 Å` — atomic-fusion only. 100% pass it; half still
   clash by vdW. **Lever:** add a van-der-Waals inter-fragment clash term to acceptance/scan.
2. **The FF clean has no whole-complex steric term.** The RDKit clean force-fields *ligand atoms
   only*, metal as fixed position constraints — two ligands can be driven together with no
   restoring force. **Lever:** a global (inter-ligand) non-bonded term.
3. **Structures are never physically relaxed.** The whole sweep is FF-only
   (`optimizer_effective:"ff"`, `xtb_available:false`). **Lever:** wire the already-coded
   `generator3d/ml_optimizer.ASEOptimizer` (xtb/MACE) as a final relax tier — needs the binary
   installed, which would also unlock a *valid* relaxation-quality metric.
4. **Conformer selection is steric-blind.** Distortion concentrates in high-CN / low-symmetry
   geometries (PBP, SQA, TBP, SPY) and big early metals (Y, Sc, Zr, Ti, Ta) — consistent with
   donor placement onto ideal template axes ignoring inter-ligand bulk. **Lever:** score
   candidate embeddings by whole-complex clash, not only coordination-sphere fit.
5. **The round-trip gate is geometry-blind.** Topological successes clash as much as failures.
   **Lever:** add a physical-quality gate (clash/strain or full-molecule RMSD) beside the
   string/RMSD gates so bad geometry can't ship as a "pass".

## Metric definitions

| Metric | Definition | Dir. | Bounds L→U |
|---|---|---|---|
| clash_vdw | # non-bonded, non-geminal pairs with dist < 0.75·(rvdw_i+rvdw_j) | min | 0→6 |
| clash_severe | same, threshold 0.60 | min | 0→2 |
| worst_overlap | min over non-bonded pairs of dist/(rvdw_i+rvdw_j) | max | 0.60→0.90 |
| angle_strain | mean \|angle − ideal(by degree)\| at non-metal centers (deg) | min | 4→20 |
| cn_divergence | \|perceived metal CN (gen) − perceived metal CN (input)\| | min | 0→2 |
| bondcount_divergence | \|perceived bond count (gen) − perceived bond count (input)\| | min | 0→6 |
| perception_ok | 1 if the sweep re-encoded the generated coords (smiles_2 present) else 0 | max | 0→1 |
| coord_rmsd | harness coordination-sphere RMSD to input (sentinels ≥900 → N/A) | min | 0.30→1.50 |
| full_divergence | symmetric element-aware heavy-atom overlay divergence to input (Å) | min | 0.50→3.00 |

Connectivity is perceived by a distance/covalent-radius cutoff (1.3·ΣR_cov); vdW/covalent radii
from RDKit; metal identity from the OIN metal tag. CN and bond-count divergence compare
generated-perceived vs input-perceived (same method both sides), so hapticity over-counting
cancels rather than becoming a false signal.

## Reproducing

```
# full non-quick capstone population (+ ligands-only UFF proxy; a few minutes)
uv run python tools/structure_distortion_report.py \
  --results-dir tmCAT-tmPHOTO_xyz_dataset/results-capstone-v042 \
  --output-dir results-v043-distortion

# fast pass without the (informational) UFF proxy, or a diverse sample
uv run python tools/structure_distortion_report.py --no-uff --output-dir <out>
uv run python tools/structure_distortion_report.py --sample 40 --output-dir <out>
```

Outputs `distortion_metrics.json` (per-molecule raw metrics + desirabilities + angle scores +
MPO; `config` block records the exact bounds/weights) and `distortion_report.md` (these tables).
Retune scoring without code edits via `--config <json>` overriding `angle_weights` / `metrics`.
