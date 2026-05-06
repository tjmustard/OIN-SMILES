# OIN-SMILES

**OIN-SMILES** is a standalone, open-source Python library for lossless conversion between 3D XYZ structures and 1D SMILES strings for Transition Metal Complexes (TMCs). It implements the **Open Isomer Notation (OIN)** to preserve stereochemical fidelity.

## Features

- **Lossless Round-Tripping**: XYZ → OIN → XYZ with exact isomer preservation.
- **Open Isomer Notation (OIN) v3.6**: Compact inline format encoding coordination geometry, slot assignments, hapticity, winding direction, and P/N stereochemistry.
- **Robust Graph Generation**: Powered by the Jensen Group's `xyz2mol` algorithm for TMCs.
- **Deterministic 3D Generation**: Uses **SCINE Molassembler** as backend — template-based placement for all ligand types, distance geometry (DG) for fallback conformer generation. Special handling for aromatic η-ligands (Cp, indenyl) via ETKDG embedding with de-aromatization to avoid RDKit kekulization failures.
- **P/N Stereocenter Encoding**: 3D-derived CIP codes for chiral phosphorus and nitrogen centers are encoded directly in the OIN string.
- **CLI**: `oin-smiles` command for one-line conversions.

## Installation

This project uses `uv` for dependency management.

1. **Install `uv`** (if not already installed):

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

2. **Clone the repository:**

    ```bash
    git clone https://github.com/tjmustard/OIN-SMILES.git
    cd OIN-SMILES
    ```

3. **Sync dependencies:**

    ```bash
    uv sync
    ```

## Usage

### CLI

```bash
# XYZ → OIN
oin-smiles xyz2oin complex.xyz

# OIN → XYZ (prints XYZ block to stdout)
oin-smiles oin2xyz "[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
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

generator = OIN3DGenerator()
result = generator.generate("[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}")

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

`generate()` accepts an optional `timeout` (seconds, default 60) passed to the `OIN3DGenerator` constructor.

### OIN v3.6 Inline Format

**Cisplatin** (square planar, *cis* isomer):

```text
[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}
```

**Transplatin** (square planar, *trans* isomer):

```text
[Pt@SP2_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}
```

**Key fields:**

| Element | Meaning |
|---|---|
| `[Pt@SP1_SPL]` | Metal atom; `SP1` = stereo descriptor (atom-order dependent); `SPL` = Square Planar geometry template |
| `{n}` | Slot tag — assigns the preceding binding atom to slot *n* of the geometry template |
| `{n>}` / `{n<}` | Slot tag with winding direction (CW / CCW) for η-ligands |
| `.` | Fragment separator (metal + each ligand) |

**Geometry templates:** `LIN`, `TPL`, `TET`, `SPL`, `SPY`, `TBP`, `OCT`, `PBP`, `TPY`

### Timeout Configuration

```python
generator = OIN3DGenerator(timeout=120)  # 120-second DG timeout
```

A `MolassemblerTimeoutError` is raised if generation exceeds the timeout.

## Development

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
uv run python tests/integration/verify_roundtrip.py [--output-dir /tmp/results]

# Compare DG strategies (single / ensemble / directed) side-by-side
uv run python tests/integration/compare_dg_strategies.py [--output-dir /tmp/results]
```

All scripts write named output artifacts when `--output-dir` is specified:
- `Ex1_CisPlatinXYZ-OIN-SMILES_original.xyz` / `.mol` / `.sdf` — input structure
- `Ex1_CisPlatinXYZ-OIN-SMILES_generated.xyz` / `.mol` / `.sdf` — generated structure
- `.mol` and `.sdf` files include full bond connectivity for generated structures

Generated MOL/SDF files use **V3000 format** when dative bonds are present (M–L bonds), which is standard for organometallic structures.

### Verified Examples

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

## Acknowledgements

- **SCINE Molassembler** — 3D structure generation and distance geometry.
- **xyz2mol** — Jensen Group's algorithm for robust graph generation from 3D coordinates.
- **OpenBabel & XTB** — Chemical file handling and geometry optimization.

## License

MIT
