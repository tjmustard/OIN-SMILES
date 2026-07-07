<div align="center">
    <img src="./media/OIN-SMILES-logo-dark.webp" alt="OIN-SMILES Logo" width="600"/>
    <h1>OIN-SMILES</h1>
    <h3><em>Lossless XYZ ↔ SMILES conversion for Transition Metal Complexes.</em></h3>
</div>

<p align="center">
    <strong>An open-source Python library implementing Open Isomer Notation (OIN) — a 1D string format that captures coordination geometry, slot assignments, hapticity, winding direction, and P/N stereochemistry for perfect 3D reconstruction.</strong>
</p>

<p align="center">
    <a href="https://github.com/tjmustard/OIN-SMILES/releases/tag/v0.3.3"><img src="https://img.shields.io/badge/release-v0.3.3-blue" alt="Latest Release"/></a>
    <a href="https://github.com/tjmustard/OIN-SMILES/actions/workflows/ci.yml"><img src="https://github.com/tjmustard/OIN-SMILES/actions/workflows/ci.yml/badge.svg" alt="CI"/></a>
    <a href="https://github.com/tjmustard/OIN-SMILES/stargazers"><img src="https://img.shields.io/github/stars/tjmustard/OIN-SMILES?style=social" alt="GitHub Stars"/></a>
    <a href="https://github.com/tjmustard/OIN-SMILES/blob/main/LICENSE"><img src="https://img.shields.io/github/license/tjmustard/OIN-SMILES" alt="License"/></a>
</p>

---

## Table of Contents

- [🔬 What is OIN-SMILES?](#-what-is-oin-smiles)
- [✨ Features](#-features)
- [⚡ Installation](#-installation)
- [🚀 Usage](#-usage)
- [✅ Verified Examples](#-verified-examples)
- [🛠️ Development](#️-development)
- [🙏 Acknowledgements](#-acknowledgements)
- [📄 License](#-license)

## 🔬 What is OIN-SMILES?

Standard SMILES notation is lossy for transition metal complexes (TMCs): coordination geometry, isomer identity, and P/N stereochemistry are destroyed when you convert a 3D XYZ structure to a 1D string. This breaks any workflow that requires exact isomer preservation — database curation, ML training sets, round-trip structure exchange.

**OIN-SMILES** solves this by implementing **Open Isomer Notation (OIN)** — an extended SMILES format that encodes everything needed to reconstruct the original 3D geometry from the string alone. A cisplatin complex becomes `[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}`; that string unambiguously reconstructs the *cis* square planar geometry, not the *trans* isomer.

## ✨ Features

- **Lossless Round-Tripping**: XYZ → OIN → XYZ → OIN with exact isomer preservation.
- **Open Isomer Notation (OIN) v3.7**: Compact inline format encoding coordination geometry, slot assignments, hapticity, winding direction, and P/N stereochemistry. The metal token is descriptor-free (`[Pt_SPL]`); cis/trans and fac/mer isomerism is carried entirely by slot order. Parsers still accept legacy `@desc` tokens.
- **Robust Graph Generation**: Powered by the Jensen Group's `xyz2mol` algorithm for TMCs.
- **3D Generation**: The vendored **MetalloGen** engine is the default backend (dummy-metal + RDKit `CoordMap` embed, constrained MMFF/UFF cleanup, standard `g-xTB` refinement with optional MACE accuracy enhancement; uses coordination-geometry-matched conformer selection); **SCINE Molassembler** remains available as the `legacy` backend (template placement + distance-geometry fallback, and the reference for Zone-A P stereo enforcement). Special handling for aromatic η-ligands (Cp, indenyl) via ETKDG embedding with de-aromatization to avoid RDKit kekulization failures.
- **P Stereocenter Encoding & Enforcement**: metal-bound chiral phosphorus centers are encoded as `[P@]`/`[P@@]` from the 3D structure and generate the correct enantiomer on both tetrahedral and square-planar complexes. (Nitrogen stereocenters are carried on backbone atoms; direct `[N@]` encoding is deferred — see CHANGELOG.)
- **Haptic-Face Round-Tripping**: η-ligand winding markers (`{n>}`/`{n<}`) survive the round trip and control which ring face the metal binds during 3D generation.
- **CLI**: `oin-smiles` command for one-line conversions.

## ⚡ Installation

This project uses [`uv`](https://docs.astral.sh/uv/) for dependency management. The
default install is lightweight — the fast FF + g-xTB path, with **no `torch`**. The
machine-learning (MACE) optimizer is an opt-in extra.

1. **Install `uv`** (if not already installed):

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

2. **Clone the repository:**

    ```bash
    git clone https://github.com/tjmustard/OIN-SMILES.git
    cd OIN-SMILES
    ```

3. **Install dependencies:**

    ```bash
    uv sync
    ```

4. **Install `g-xTB`** (recommended — the default 3D-generation optimizer):

    OIN-SMILES defaults to Grimme Lab's `g-xTB` for fast, accurate structural refinement.
    ```bash
    bash tools/install_gxtb.sh
    ```
    *Detects your OS/architecture (Linux/macOS) and drops the static binary into
    `.venv/bin` to keep your system clean. If `xtb` is absent, 3D generation falls
    back to the force field automatically.*

5. **(Optional) MACE machine-learning optimizer** — highest accuracy:

    ```bash
    uv sync --extra mace                  # pulls mace-torch + a pinned CUDA 11.8 torch
    bash tools/install_mace_weights.sh    # downloads the weights + writes .env
    ```
    See [`docs/OPTIMIZERS.md`](docs/OPTIMIZERS.md) and
    [`models/mace/README.md`](models/mace/README.md) for details.

> [!NOTE]
> The `mace` extra pulls PyTorch with CUDA 11.8 (configured via `tool.uv.index` in
> `pyproject.toml`). For a CPU-only or different-CUDA build, adjust that index URL
> before syncing. The default `uv sync` installs neither `torch` nor `mace-torch`.

## 🚀 Usage

### CLI

```bash
# XYZ → OIN
uv run oin-smiles xyz2oin complex.xyz

# OIN → XYZ (prints XYZ block to stdout).
# Default backend is the MetalloGen engine refined with standard g-xTB.
uv run oin-smiles oin2xyz "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"

# Higher accuracy MACE refinement (needs the `mace` extra + weights; note --extra mace):
uv run --extra mace oin-smiles oin2xyz "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}" --optimizer mace-omol-0-extra-large-1024

# Fast FF-only path (no torch/xtb), or the legacy Molassembler backend:
uv run oin-smiles oin2xyz "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}" --optimizer ff
uv run oin-smiles oin2xyz "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}" --engine legacy
```

### 3D to 1D (XYZ to OIN)

```python
from oinsmiles import XYZToSMILES

oin = XYZToSMILES().convert("complex.xyz")
print(oin)
```

### 1D to 3D (OIN to XYZ)

```python
from oinsmiles.generation.engine import OIN3DGenerator

# Default engine is "metallogen" refined with standard g-xTB (requires g-xTB binary).
# Use optimizer="mace-omol-0-extra-large-1024" for higher accuracy MLIP refinement.
# Use optimizer="ff" for the fast FF-only path, or engine="legacy"
# for the Molassembler backend (the reference for Zone-A P stereo enforcement).
generator = OIN3DGenerator()
result = generator.generate("[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}")

# XYZ block (always available)
with open("generated.xyz", "w") as f:
    f.write(result.xyz)

# Bonded RDKit mol — includes M–L dative bonds and full ligand connectivity
# (None if bond connectivity could not be determined, e.g. for eta fallbacks)
if result.mol is not None:
    from rdkit import Chem
    Chem.MolToMolFile(result.mol, "generated.mol")
    writer = Chem.SDWriter("generated.sdf")
    writer.write(result.mol)
    writer.close()
```

### 3D Optimization Workflow

The MetalloGen engine implements a multi-stage optimization pipeline for 1D → 3D generation:

1. **Force Field Relaxation (Default)**: Uses constrained MMFF/UFF to relax the generated conformers into the geometric template. Convergence criteria are controlled by `ff_preset` (`loose`, `default`, `tight`, `very_tight`).

2. **MLIP / Semi-empirical Refinement (Optional)**: Refines the FF-relaxed geometry pool using an advanced optimizer like standard `g-xTB` (default) or MACE (e.g., `mace-omol-0-extra-large-1024`), then energy-ranks the ensemble. The returned conformer is chosen by **coordination-geometry match** to the requested template (e.g. square-planar), with energy breaking ties — so a floppy donor that admits an energetically competitive distorted geometry still yields the correct isomer. Haptic (η) ligands fall back to lowest-energy.

```python
# Generate an ensemble of 5, pre-optimize with tight FF, and refine with MACE
generator = OIN3DGenerator(
    ff_preset="tight",
    optimizer="mace-omol-0-extra-large-1024",
    ensemble_size=5,
    timeout=120
)
```

A `MolassemblerTimeoutError` is raised if generation exceeds the `timeout` limit (seconds).

### OIN v3.7 Inline Format

**Cisplatin** (square planar, *cis* isomer):

```text
[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}
```

**Transplatin** (square planar, *trans* isomer):

```text
[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}
```

**Key fields:**

| Element | Meaning |
|---|---|
| `[Pt_SPL]` | Metal atom + geometry template (`SPL` = Square Planar); isomerism (cis/trans, fac/mer) is carried entirely by slot order, not by the metal token |
| `{n}` | Slot tag — assigns the preceding binding atom to slot *n* of the geometry template |
| `{n>}` / `{n<}` | Slot tag with winding direction (CW / CCW) for η-ligands |
| `[P@]` / `[P@@]` | Metal-bound phosphorus stereocenter (lone-pair CIP convention) |
| `.` | Fragment separator (metal + each ligand) |

**Geometry templates:** `LIN`, `TPL`, `TET`, `SPL`, `SPY`, `TBP`, `OCT`, `PBP`, `TPY`



## ✅ Verified Examples

The following complexes pass the full round-trip test (OIN string identity + RMSD < 1.0 Å):

**Square Planar & Simple Geometries**
- **Cisplatin** — square planar Pt
- **Transplatin** — square planar Pt, trans isomer
- **Ferrocene** — metallocene with eclipsed Cp rings
- **VOacac2** — square pyramidal vanadium with bidentate acetylacetonate ligands

**P/N Stereocenter Encoding**
- **PdCl2-R-BINAP** — square planar Pd with axial-chiral BINAP diphosphine (2,2'-bis(diphenylphosphino)-1,1'-binaphthyl)
- **PdCl2-RR-BDNN** — square planar Pd with RR-BDNN N-chiral diphosphine (N,N-bis(diethylamino)naphthalene with chiral centers on nitrogen)
- **PdCl2-RR-BDPP** — square planar Pd with RR-BDPP P-chiral diphosphine (phosphorus stereocenters with proper @/@@ encoding in SMILES)

**Ansa-Metallocenes (Aromatic η-Ligands)**
- **TiCat1** — tetrahedral Ti with bridged Cp–Si(Me)₂–Cp ligand (3D generation fixed; round-trip bonding inference incomplete due to de-aromatization requirement)
- **TiCat3** — tetrahedral Ti with bridged indenyl–Si(Me)₂–indenyl ligand (3D generation fixed)
- **TiCat4** — tetrahedral Ti with bridged indenyl–Si(Me)₂–indenyl ligand variant (3D generation fixed)

## 🛠️ Development

### Running Tests

```bash
uv run python -m unittest discover tests
```

### Verification Scripts

Three scripts are available for QA:

```bash
# Fast smoke test (~4 examples, no tmQM dataset)
bash tests/run_verification_fast.sh

# Full suite — all manual examples
bash tests/run_verification.sh

# Comprehensive — all examples including tmQM dataset (~103 complexes)
bash tests/run_verification_ALL.sh
```

Or run individual scripts:

```bash
# XYZ → OIN correctness
uv run python tests/integration/verify_xyz_to_oin.py [--include-tmqm]

# Full round-trip: XYZ → OIN → XYZ → OIN (RMSD + string identity)
# Writes XYZ, MOL, SDF and OIN files to --output-dir when provided
uv run python tests/integration/verify_roundtrip.py [--output-dir /tmp/results] [--optimizer <opt>] [--ff-preset <preset>]

# Compare DG strategies (single / ensemble / directed) side-by-side
uv run python tests/integration/compare_dg_strategies.py [--output-dir /tmp/results]
```

All scripts write named output artifacts when `--output-dir` is specified:
- `Ex1_CisPlatinXYZ-OIN-SMILES_original.xyz` / `.mol` / `.sdf` — input structure
- `Ex1_CisPlatinXYZ-OIN-SMILES_generated.xyz` / `.mol` / `.sdf` — generated structure
- `.mol` and `.sdf` files include full bond connectivity for generated structures

> [!NOTE]
> Generated MOL/SDF files use **V3000 format** when dative bonds are present (M–L bonds), which is standard for organometallic structures.

## 🙏 Acknowledgements

- **[SCINE Molassembler](https://github.com/qcscine/molassembler/)** — 3D structure generation and distance geometry.
- **[MetalloGen](https://github.com/kyunghoonlee777/MetalloGen)** — 3D generation engine and optimization workflows for transition metal complexes.
- **[MACE](https://github.com/acesuit/mace)** — Fast and accurate Machine Learning Interatomic Potentials (MLIP) for 3D geometry refinement.
- **[xyz2mol](https://github.com/jensengroup/xyz2mol_tm)** — Jensen Group's algorithm for robust graph generation from 3D coordinates.
- **[OpenBabel](https://github.com/openbabel/openbabel) & [XTB](https://github.com/grimme-lab/xtb)** — Chemical file handling and semi-empirical geometry optimization.

## 📄 License

This project is licensed under the terms of the MIT open source license. Please refer to the [LICENSE](./LICENSE) file for the full terms.
