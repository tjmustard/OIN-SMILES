# tmCAT / tmPHOTO — Functional Transition Metal Complex Dataset

> **A homage to Heather J. Kulik and the Kulik Group (MIT).**
> This branch (`Kulik_TMC_Dataset`) is a data-only distribution of the **tmCAT** and
> **tmPHOTO** subsets of functional transition metal complexes, curated by the Kulik group
> and used here as the round-trip validation corpus for the **OIN-SMILES** project. It is
> shared as a shout-out to their work generating this dataset.

---

## Credit & Citation

These structures originate from the datasets described in:

> Kevlishvili, I.; St. Michel, R. G.; Garrison, A. G.; Toney, J. W.; Adamji, H.; Jia, H.;
> Román-Leshkov, Y.; **Kulik, H. J.** *"Leveraging natural language processing to curate the
> tmCAT, tmPHOTO, tmBIO, and tmSCO datasets of functional transition metal complexes."*
> **Faraday Discussions**, **2025**, *256*, 275–303. DOI: [10.1039/d4fd00087k](https://doi.org/10.1039/d4fd00087k)

The source work defines **four** subsets of functional transition metal complexes
(**tmCAT** = catalysis, **tmPHOTO** = photochemistry, **tmBIO** = biological, **tmSCO** =
spin-crossover), curated from the literature with natural-language-processing methods and
built on Cambridge Structural Database (CSD) entries. **This branch ships only the two
subsets used by OIN-SMILES: `tmCAT` and `tmPHOTO`.**

If you use these structures, please cite the paper above and acknowledge the Kulik Group.
Underlying crystal structures are subject to the CSD's terms; consult the source publication
for licensing of the derived dataset.

---

## What's in this branch

```
tmCAT-tmPHOTO_xyz_dataset/
├── README.md            ← this file
├── rebuild.sh           ← helper: rebuild summary_roundtrip.json from individual reports
├── run_continuous.sh    ← helper: continuously round-trip random molecules
├── cat/    *.xyz        ← tmCAT subset  — 21,631 structures
└── photo/  *.xyz        ← tmPHOTO subset —  4,599 structures
```

- **26,230** total `.xyz` files (**25,197** unique `<REFCODE>_comp_<N>` entries; **1,033**
  refcodes appear in *both* folders because a complex can belong to more than one subset).
- **`cat/`** and **`photo/`** are flat folders — one `.xyz` per file, no per-molecule
  subdirectories.

> **This is a data-only branch.** It intentionally does **not** contain the `oinsmiles`
> library or the test harness. The scripts here are kept for reference/provenance and are
> meant to be run from a full checkout of `main` (see "Running the round-trip test" below).

---

## File format

Each file is a standard XYZ geometry with a descriptive comment line:

```
85
Refcode: ABAFOZ_comp_0 | Dataset: tmCAT | Charge: 0 | Potential Spins: [1, 3, 5]
Pd      3.145400     4.235200     3.798300
...
```

- **Line 1** — atom count.
- **Line 2** — metadata: `Refcode`, `Dataset`, `Charge`, `Potential Spins`.
  - `Dataset` may list one or more subsets, e.g. `tmCAT`, `tmPHOTO`, or `tmCAT;tmPHOTO`
    (semicolon-separated) when the complex belongs to both. **The containing folder
    (`cat/` vs `photo/`) is the authoritative subset assignment for this branch.**
  - `Charge` is the total molecular charge; `Potential Spins` lists candidate spin
    multiplicities.
- **Lines 3+** — `element  x  y  z` (Ångström).

**Naming:** `<REFCODE>_comp_<N>.xyz`, where `<REFCODE>` is the 6-letter CSD refcode and
`<N>` is the component/conformer index (the vast majority are `_comp_0`; a minority have
`_comp_1` … `_comp_12`).

---

## How OIN-SMILES uses this dataset

OIN-SMILES performs **lossless conversion between 3D XYZ structures and 1D OIN-SMILES
strings** for transition metal complexes. This corpus is the large-scale benchmark for that
conversion: every structure is used to validate the **XYZ → OIN → 3D → OIN** round-trip,
checking that the canonical notation, coordination geometry (RMSD), and atom counts survive
the trip. Because these are real, chemically diverse TMCs (varied metals, coordination
numbers, chelates, η-bonded ligands, macrocycles), they stress-test the encoder/generator far
beyond hand-built fixtures.

---

## Running the round-trip test

The harness lives on the project's **`main`** branch, not here. From a full checkout of the
repository (with the `oinsmiles` package installed via `uv sync`), point the harness at this
dataset directory:

```bash
# from the repository root on `main`, with this dataset present as tmCAT-tmPHOTO_xyz_dataset/
uv run python tools/test_dataset_roundtrip.py \
    --dataset-dir tmCAT-tmPHOTO_xyz_dataset \
    --output-dir  <some-results-dir> \
    --quick            # 60s/mol fast pass; drop for the full g-xTB pass
```

Useful flags: `--limit N`, `--random`, `--continue` (skip already-processed molecules),
`--rerun-failed`, `--only "A,B,C"`, `--shard I:N` (parallel workers), `--mol-timeout SECONDS`
(kill hung molecules). The harness walks the dataset directory for every `*.xyz` (ignoring
`*_generated.xyz` and anything inside the output dir), so `cat/` and `photo/` are picked up
automatically.

**Included helper scripts** (run from the repo root on a full checkout; they `cd ..` and call
`uv run`):

- `run_continuous.sh` — infinite loop that round-trips one random, not-yet-processed molecule
  at a time (`--quick --limit 1 --random --continue`).
- `rebuild.sh` — rebuilds `summary_roundtrip.json` from the per-molecule reports
  (`tools/rebuild_summary.py`).

Results (`summary_roundtrip.json`, `individual_reports/`, `structures/`, `test_failures/`)
are written to the chosen output directory, which is kept **out of version control** — this
branch ships only the input structures, not run output.

---

## For agents

- This branch contains **only** the dataset folder and two helper scripts — no `src/`,
  `tools/`, or `tests/`. Do not expect to import `oinsmiles` from here.
- `cat/` and `photo/` are flat directories of `.xyz` files; molecule identity is the file
  basename (`<REFCODE>_comp_<N>`).
- To actually run the round-trip harness, switch to `main` (which has the library +
  `tools/test_dataset_roundtrip.py`) and pass `--dataset-dir` pointing at this folder.
- The prior run outputs (`20260707-results/`, `results-archive-20260706/`) and the
  scratch script `run_continuous-v0.4.0.sh` were **deliberately excluded** from this branch.
