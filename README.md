# OIN-SMILES

**OIN-SMILES** is a standalone, open-source Python library capable of lossless conversion between 3D XYZ structures and 1D SMILES strings for Transition Metal Complexes (TMCs). It implements the **Open Isomer Notation (OIN)** to preserve stereochemical fidelity.

## Features

- **Lossless Round-Tripping**: XYZ -> SMILES -> XYZ with exact isomer preservation using OIN.
- **Open Isomer Notation (OIN) V3.6**: Implements the V3.6 Inline Syntax specifically designed for reliability and readability.
- **Explicit Directionality**: Introduces directional winding tags (`<`, `>`) to correctly handle haptic ligand faces (e.g., distinguishing between top-face and bottom-face coordination for Cp rings).
- **Robust Graph Generation**: Powered by the Jensen Group's `xyz2mol` algorithm for Transition Metal Complexes.
- **Deterministic 3D Generation**: Uses **Architector** as a backend to reconstruct 3D structures from OIN strings, enforcing specific coordination geometries defined by the OIN vectors.

## Installation

This project relies on several scientific packages that are best managed via Conda. Note that this environment installs a specific development branch of **Architector** directly from GitHub to support secondary solvation shells.

### 1. Prerequisite: Install Micromamba
If you do not have a Conda-based manager, install Micromamba:
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

### 2. Create the Environment
Clone the repository and build the environment using the provided YAML file.

git clone https://github.com/tjmustard/OIN-SMILES.git
cd OIN-SMILES

# Create environment (installs xtb, rdkit, and the custom Architector branch)
micromamba create -f environment.yml

# Activate
micromamba activate oin-smiles

### 3. Verify
python -c "from ase.calculators.xtb import XTB; import rdkit; import architector; print('Environment is ready.')"

## Directory Structure

```text
.
├── PRDs/                           # Product Requirements Documents
├── src/
│   └── oinsmiles/
│       ├── core/                   # Core functionality (Graph, Translator)
│       ├── generation/             # 1D to 3D Generation Logic
│       │   ├── architector_adapter.py # Adapter for Architector
│       │   ├── engine.py           # Generation engine entry point
│       │   └── oin_parser.py       # Parser for OIN generation
│       ├── oin/                    # OIN Syntax Handling
│       │   ├── inline.py           # Inline syntax V3.6 processing
│       │   ├── parser.py           # OIN string parsing
│       │   └── writer.py           # OIN string writing
│       └── utils/                  # Utilities (xyz2mol, aligners)
├── tests/                          # Integration and Unit Tests
├── LICENSE
├── pyproject.toml
└── README.md
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

### Example OIN-SMILES (V3.6)

**Cisplatin:**

```text
[Pt_SPL].[Cl]{0}.[Cl]{1}.[NH3]{2}.[NH3]{3}
```

**Ferrocene (Eclipsed):**

```text
[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1
```

**Syntax Guide:**

1.  **Tag Injection**: Tags are now injected directly into the SMILES string (Inline Syntax) rather than in a separate block.
    *   **Metal Geometry**: `[Pt_SPL]` (Metal Symbol + Geometry Code).
    *   **Ligand Slots**: `{0}`, `{1}`, etc. attached to the binding atom(s).
2.  **Explicit Directionality**:
    *   `{n>}`: Indicates **Forward/Clockwise** winding relative to the template face.
    *   `{n<}`: Indicates **Reverse/Counter-Clockwise** winding.
    *   This is critical for distinguishing haptic faces (e.g. top vs bottom of a Cp ring).
3.  **Sanitization**: Zone A atoms (bonded to metal) are locked with explicit Hydrogen counts (e.g. `[NH3]`, `[cH]`) to prevent "SMILES Drift".
4.  **Strict Bracketing**: Coordinating atoms (C, O, H) are explicitly bracketed (e.g., `[C]#O`, `[H]`) to ensure correct graph reconstruction and avoid valence ambiguity.

### 1D to 3D (OIN to XYZ)

Generate a 3D structure from an OIN string. This uses the `OIN3DGenerator` which leverages `Architector` to build the complex.

```python
from oinsmiles.generation.engine import OIN3DGenerator

generator = OIN3DGenerator()
# V3.6 String for Cisplatin
oin_string = "[Pt_SPL].[Cl]{0}.[Cl]{1}.[NH3]{2}.[NH3]{3}"

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

### Real-Life Examples

To run a set of real-life examples:

```bash
uv run python tests/integration/verify_xyz_to_oin.py
```

#### Verified Examples (V3.6)

The following complex examples have been verified to work correctly with OIN-SMILES V3.6:

- **Cisplatin**: Square planar platinum complex.
- **Transplatin**: Trans isomer of Cisplatin.
- **Ferrocene**: Metallocene with eclipsed Cp rings (Validates Haptic Directionality).
- **fac-Ir(ppy)3**: Facial isomer of Iridium tris(phenylpyridine).
- **mer-Ir(ppy)3**: Meridional isomer.
- **TiCp2Me2**: Titanocene dimethyl (Validates Mixed Haptic/Monodentate).
- **TiCat Series**: Complex catalysts demonstrating advanced haptic winding.
- **Hydrides & Carbonyls**: Correctly brackets binding atoms (e.g., `[H]`, `[C]#O`) to prevent valence errors.

## Acknowledgements

- **Architector**: Used for 3D structure generation and assembly.
- **xyz2mol**: The Jensen Group's algorithm is used for robust graph generation from 3D coordinates.
- **OpenBabel & XTB**: Used for chemical file handling and geometry optimization.

## License

MIT
