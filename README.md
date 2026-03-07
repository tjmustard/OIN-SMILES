# OIN-SMILES

**OIN-SMILES** is a standalone, open-source Python library for lossless conversion between 3D XYZ structures and 1D SMILES strings for Transition Metal Complexes (TMCs). It implements the **Open Isomer Notation (OIN)** to preserve stereochemical fidelity.

## Features

- **Lossless Round-Tripping**: XYZ → OIN → XYZ with exact isomer preservation.
- **Open Isomer Notation (OIN) v3.6**: Compact inline format encoding coordination geometry, slot assignments, hapticity, winding direction, and P/N stereochemistry.
- **Robust Graph Generation**: Powered by the Jensen Group's `xyz2mol` algorithm for TMCs.
- **Deterministic 3D Generation**: Uses **SCINE Molassembler** as backend — template-based placement for all ligand types, distance geometry (DG) for fallback conformer generation.
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
xyz_block = generator.generate("[Pt@SP1_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}")

with open("generated.xyz", "w") as f:
    f.write(xyz_block)
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

```bash
# XYZ → OIN correctness
uv run python tests/integration/verify_xyz_to_oin.py

# Full round-trip: XYZ → OIN → XYZ → OIN (RMSD + string identity)
uv run python tests/integration/verify_roundtrip.py
```

### Verified Examples

The following complexes pass the full round-trip test (OIN string identity + RMSD < 1.0 Å):

- **Cisplatin** — square planar Pt
- **Transplatin** — square planar Pt, trans isomer
- **Ferrocene** — metallocene with eclipsed Cp rings
- **VOacac2** — square pyramidal vanadium with bidentate acetylacetonate ligands

## Acknowledgements

- **SCINE Molassembler** — 3D structure generation and distance geometry.
- **xyz2mol** — Jensen Group's algorithm for robust graph generation from 3D coordinates.
- **OpenBabel & XTB** — Chemical file handling and geometry optimization.

## License

MIT
