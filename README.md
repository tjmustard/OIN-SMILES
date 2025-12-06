# OIN-SMILES

**OIN-SMILES** is a standalone, open-source Python library capable of lossless conversion between 3D XYZ structures and 1D SMILES strings for Transition Metal Complexes (TMCs). It implements the **Open Isomer Notation (OIN)** to preserve stereochemical fidelity.

## Features

- **Lossless Round-Tripping**: XYZ -> SMILES -> XYZ with exact isomer preservation using OIN.
- **Open Isomer Notation (OIN)**: Human-readable metadata tags for explicit relative coordinates, dative bonds, and hapticity.
- **Robust Graph Generation**: Powered by the Jensen Group's `xyz2mol` algorithm for Transition Metal Complexes.
- **Deterministic 3D Generation**: Uses **Architector** as a backend to reconstruct 3D structures from OIN strings, enforcing specific coordination geometries defined by the OIN vectors.

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

### 3D to 1D (XYZ to SMILES)

Convert an XYZ file into an OIN-SMILES string that captures the exact stereochemistry.

```python
from oinsmiles import XYZToSMILES

converter = XYZToSMILES()
smiles = converter.convert("complex.xyz")
# Returns an OIN-SMILES string containing connectivity and 3D metadata
print(smiles)
```

### Example OIN-SMILES (V2.4)

**Cisplatin:**

```text
[Pt].[Cl].[Cl].[NH3].[NH3] |g:SPL|w:1.0:0;2.0:1;3.0:2;4.0:3|
```

**Key Features:**

-   **Compact SMILES**: Uses "Sanitization-First" to generate canonical, compact SMILES strings (e.g., `[NH3]` instead of `[H]N([H])[H]`). Explicit Hydrogen counts are enforced on Zone A atoms (bonded to metal).
-   **Disconnected SMILES**: Metal and ligands are separated by dots `.`
-   **Geometry Tag (g)**: Explicitly defines the coordination geometry template (e.g., `g:SPL` for Square Planar).
-   **Index-Based Topology (w)**: Defines the geometry using discrete slot assignments.
    -   Format: `w:Rank.Idx:Slot;...`
    -   `Rank`: 1-based Ligand Index in the dot-separated SMILES string (Metal is 0).
    -   `Idx`: 0-based Atom Index within the ligand's canonical SMILES string.
    -   `Slot`: 0-based Slot Index on the Geometry Template.
-   **Deterministic**: Maximizes vector alignment to standard templates for a unique canonical string.

### 1D to 3D (OIN to XYZ)

Generate a 3D structure from an OIN string. This uses the `OIN3DGenerator` which leverages `Architector` to build the complex.

```python
from oinsmiles.generation.engine import OIN3DGenerator

generator = OIN3DGenerator()
# V2.4 String for Cisplatin
oin_string = "[Pt].[Cl].[Cl].[NH3].[NH3] |g:SPL|w:1.0:0;2.0:1;3.0:2;4.0:3|"

# Generate the structure (returns an Architector Molecule object)
structure = generator.generate(oin_string)

# Write to XYZ file
structure.write_file("generated_cisplatin.xyz")
```

## Development

### Running Tests

To run the unit tests:

```bash
uv run python -m unittest discover tests
```

### Verification Scripts

We provide scripts to verify the end-to-end functionality, including round-trip tests (XYZ -> OIN -> XYZ).

```bash
# Run the full verification suite
bash tests/run_verification.sh
```

This script runs:

1. **Phase 1 Verification**: Checks basic generation capabilities.
2. **Round-Trip Verification**:
    - **Test A**: Starts with an XYZ file, converts to OIN, regenerates XYZ, and compares structures (RMSD).
    - **Test B**: Starts with an OIN string, generates XYZ, converts back to OIN, and compares the strings.

### Real-Life Examples

To run a set of real-life examples:

```bash
uv run python tests/integration/verify_xyz_to_oin.py
```

#### Verified Examples

The following complex examples have been verified to work correctly with OIN-SMILES, including accurate 3D reconstruction from the OIN string:

- **Cisplatin**: Square planar platinum complex.
- **Transplatin**: Trans isomer of Cisplatin.
- **Ferrocene**: Metallocene with eclipsed Cp rings.
- **fac-Ir(ppy)3**: Facial isomer of Iridium tris(phenylpyridine).
- **mer-Ir(ppy)3**: Meridional isomer of Iridium tris(phenylpyridine).
- **PdCl2(Butene)**: Palladium complex with an alkene ligand.
- **PdCl2(PhenPhosMe)**: Palladium complex with a phosphine ligand.

## Acknowledgements

- **Architector**: Used for 3D structure generation and assembly.
- **xyz2mol**: The Jensen Group's algorithm is used for robust graph generation from 3D coordinates.
- **OpenBabel & XTB**: Used for chemical file handling and geometry optimization.

## License

MIT
