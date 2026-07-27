# v0.4.3 A3 — vdW acceptance term (gated OFF by default)

**TL;DR.** A3 adds a whole-complex van-der-Waals steric-clash term (`generator3d/clash.py`,
a faithful lift of the release quality metric so the gate and the metric agree) at four points
that previously judged only atomic *fusion*: the embed acceptance gate, the UFF pool ranking,
the geometry selection, and the FF-clean scan. On a worst-cohort corpus sample the term drives
the generated **vdW-clash fraction from 62.5% to 2.5%** — but it does so by **ejecting bulky /
weakly-bound ligands off the metal**, which regresses the round trip. Because a geometry win
that drops round-trip pass is not a win, **the term ships gated OFF by default**
(`clash.VDW_ACCEPTANCE_ENABLED = False`); it is the bar A4's tight Kabsch placement is judged
against, and A5 decides whether to flip it on once the pool contains tight-*and*-clean
conformers.

## Why it is off by default (the finding)

The distortion research set A3 the target of the 53%-generated-vs-5%-real vdW-clash gap. The
naive fix — reject / down-rank clashing conformers at acceptance and selection — works on the
metric but has a fatal mechanism on the **current (pre-A4) conformer pool**:

> The pool contains *tight-but-clashing* conformers and *loose-but-clean* conformers, but not
> *tight-and-clean* ones. So any clash-minimizing acceptance or selection necessarily prefers
> the loose option, in which a bulky or weakly-bound ligand has splayed away from the metal.
> That geometry re-perceives as a **detached ligand**, so the round-trip OIN loses a donor.

This is real ligand ejection, not a canonical-key artifact. Examples (input → A3-on re-encode):

| molecule | input | A3-on generated re-encode |
|---|---|---|
| ABIRIO | `[Ru_TBP]` + η-diene (5-coord) | `[Ru_TPL]` — η-diene fallen off (3-coord) |
| AJEROW | `[Ru_SPY]` (5-coord) | `[Ru_TPL]` — imine-N + one P detached |
| ABOPOY | `[Zn_TPL]` (3-coord) | `[Zn_LIN]` — phosphinimine-N detached |

Real crystals are tight *and* clean, so low clash with intact coordination is achievable — but
by **better placement** (A4 Kabsch) and **relaxation** (A5 ASE), not by selecting among a pool
whose only clean members are the detached ones.

## A/B evidence

Fixed sample of 40 molecules from `results-capstone-v042` (30 worst-cohort — Zr/Ti/Y/Sc/Hf
metals and TBP/SPY/PBP geometries — plus 10 controls), regenerated pre (main `456b906a`, which
includes A0/A1/B1) vs post (main + A3 lever subsets), `optimizer="ff"`, `seed=42`. Round-trip =
`winding_canonical_key(normalize_oin_for_comparison(reencode))` matches the input OIN's key.
Scored with the release metric (`tools/structure_distortion_report.py`).

| configuration | clash fraction (≥1 vdW clash) | mean clashes | round-trip match |
|---|---|---|---|
| main (pre-A3) | 62.5% (25/40) | 2.58 | **70.0%** (28/40) |
| A3 all four levers | **2.5%** (1/40) | 0.03 | 42.5% (17/40) |
| A3 minus clean_geometry guard | 15.0% (6/40) | 0.42 | 60.0% (24/40) |
| A3 selection tiebreak only | 27.5% (11/40) | 0.75 | 62.5% (25/40) |

Every configuration that lowers clash also lowers round-trip: the two are coupled through
coordination tightness on this pool. Even the selection-only variant — restricted to conformers
that classify as the target geometry — regresses, because the least-clashing member of that set
has its donors at the edge of coordination perception.

**With the gate OFF (the shipped default) the generated structure is byte-identical to pre-A3**
(each lever's disabled branch is exactly the prior code — an import plus a flag-guarded no-op),
so the goldens keep their canonical keys and there is no round-trip regression. Verified two
ways: (1) code review — every hunk is either an import or wrapped in
`if clash.VDW_ACCEPTANCE_ENABLED:`; (2) regenerating clashy molecules flag-off and diffing
against main — byte-identical for every *deterministic* molecule (one worst-cohort SPY case,
AVAVEY, differs, but is independently non-deterministic run-to-run — same flag-off code gives
two different structures across three runs — consistent with the release's 98.2% determinism).
The full unit suite is green flag-off (`468 OK / 3 skip`).

## What ships

`generator3d/clash.py` — one shared clash definition (`vdw_clash_count`, `mol_clash_count`),
matched to the metric by a cross-check unit test. Four gated levers:

1. `embed._finalize_positions` — reject a vdW-clashing conformer and score it into `(-1, 0)` so
   the existing best-rejected fallback surfaces the least-clashing candidate.
2. `generator3d/__init__.py` — stable clash re-rank of the deduplicated UFF pool.
3. `metallogen_adapter._select_by_geometry` — clash-primary sort key among conformers that
   classify as the target geometry.
4. `clean_geometry.ff_clean` — a monotone guard that reverts an FF step which adds clash.

## Enabling it (A4 / A5)

Per run: `OIN_VDW_ACCEPTANCE=1`. In code / tests: `clash.VDW_ACCEPTANCE_ENABLED = True`. A5's
decision session should re-run the A/B above **on A4's Kabsch-placed pool** and flip the default
only if the clash fraction drops without a round-trip regression.
