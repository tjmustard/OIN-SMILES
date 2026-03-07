# OIN-SMILES

**OIN-SMILES** (v0.2.0) is a Python library for lossless conversion between 3D XYZ structures and 1D SMILES strings for Transition Metal Complexes (TMCs). It implements the **Open Isomer Notation (OIN)** to preserve stereochemical fidelity, including P/N stereocenters and axial chirality.

## Features

- **Lossless Round-Tripping**: XYZ → OIN-SMILES with exact isomer preservation.
- **OIN Inline Syntax (V3.6)**: Compact notation embedding metal geometry, ligand slots, and winding direction inside SMILES.
- **P/N Stereocenter Encoding**: Pre-fragmentation CIP assignment (`CIPAssigner`) with lone-pair pseudo-atom fallback for 3-coordinate centers.
- **Axial Chirality (`@b` tag)**: Biaryl atropisomer detection via dihedral angle encoded as a custom `@b` sidecar.
- **Molassembler Backend**: SCINE Molassembler for 3D structure generation, with `ProcessPoolExecutor` isolation and 60s timeout.
- **CLI**: `oin-smiles xyz2oin <path>` and `oin-smiles oin2xyz <oin>` for pipeline integration.
- **Robust Graph Generation**: Powered by the Jensen Group's `xyz2mol` algorithm adapted for TMCs.

## Installation

This project uses [`uv`](https://github.com/astral-sh/uv) for package management.

```bash
git clone https://github.com/tjmustard/OIN-SMILES.git
cd OIN-SMILES
uv sync
```

> **Note:** `scine-molassembler` requires the SCINE conda-forge package. See [SCINE documentation](https://scine.ethz.ch/download/molassembler) for build instructions. `xtb` bindings require an independently installed `xtb` binary.

## Directory Structure

```text
.
├── spec/
│   └── compiled/          # MiniPRDs and architecture.yml hypergraph
├── src/
│   └── oinsmiles/
│       ├── cli.py           # CLI entry point (oin-smiles)
│       ├── core/            # Translation, chirality, graph
│       │   ├── chirality.py    # CIPAssigner, ChiralityRecoveryUtility, PseudoAtomStrategy
│       │   └── translator.py   # XYZToSMILES, SMILESToXYZ
│       ├── generation/      # 3D generation
│       │   ├── engine.py          # OIN3DGenerator
│       │   ├── molassembler_adapter.py  # MolassemblerAdapter
│       │   ├── canonicalizer.py   # PAI alignment + OIN sidecar generation
│       │   └── oin_parser.py      # ParsedOIN dataclass + parser
│       ├── oin/             # OIN string parsing and writing
│       │   ├── inline.py       # V3.6 inline syntax, @b tag
│       │   ├── parser.py       # OINParser (tags, @b, bond stereo)
│       │   └── sanitizer.py    # Zone A sanitization + BondStereo recovery
│       └── utils/           # xyz2mol, aligners, sorting
├── tests/
│   ├── integration/         # .xyz fixtures + CLI tests
│   ├── unit/                # Unit tests per module
│   ├── test_regression_stability.py
│   ├── test_binap_stability.py
│   ├── test_chiral_p.py
│   ├── test_chiral_n.py
│   └── test_axial_chiral.py
├── CHANGELOG.md
├── pyproject.toml
└── README.md
```

## Usage

### XYZ → OIN-SMILES (Python)

```python
from oinsmiles import XYZToSMILES

converter = XYZToSMILES()
oin = converter.convert("complex.xyz")
print(oin)
```

### OIN-SMILES → XYZ (Python)

```python
from oinsmiles.generation.engine import OIN3DGenerator

generator = OIN3DGenerator()
xyz = generator.generate("[Pt_SPL].[Cl]{0}.[Cl]{1}.[NH3]{2}.[NH3]{3}")
print(xyz)
```

### CLI

```bash
# XYZ to OIN-SMILES
oin-smiles xyz2oin tests/integration/cisplatin.xyz

# OIN-SMILES to XYZ
oin-smiles oin2xyz "[Pt_SPL].[Cl]{0}.[Cl]{1}.[NH3]{2}.[NH3]{3}"

# Help
oin-smiles --help
```

Exit codes: `0` success, `1` generic error, `2` Molassembler timeout.

## OIN Syntax (V3.6)

| Element | Example | Meaning |
|---|---|---|
| Metal + geometry | `[Pt_SPL]` | Square planar Pt |
| Ligand slot | `{0}` | Binds to slot 0 |
| Haptic winding | `{0>}` | Forward/CW face |
| Chirality | `[P@@H]` | P stereocenter |
| Axial chirality | `@b 5-6:STEREOATROP_CCW` | Biaryl atropisomer |

**Cisplatin:**
```
[Pt_SPL].[Cl]{0}.[Cl]{1}.[NH3]{2}.[NH3]{3}
```

**Ferrocene (eclipsed):**
```
[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1
```

## Public API

```python
from oinsmiles import (
    XYZToSMILES,              # XYZ file → OIN string
    SMILESToXYZ,              # OIN string → TMCGraph
    CIPAssigner,              # Pre-fragmentation CIP assignment
    ChiralityRecoveryUtility, # Re-applies @/@@ from _CIPCode
    PseudoAtomStrategy,       # Lone-pair pseudo-atom for 3-coord P/N
)
from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.generation.molassembler_adapter import MolassemblerTimeoutError
```

## Running Tests

```bash
# All tests
PYTHONPATH=src uv run python -m unittest discover tests

# Specific suites
PYTHONPATH=src uv run python -m unittest tests/test_regression_stability.py
PYTHONPATH=src uv run python -m unittest tests/test_chiral_p.py
PYTHONPATH=src uv run python -m unittest tests/integration/test_cli.py
```

### Verified Complexes (v0.2.0)

| Complex | Geometry | Notes |
|---|---|---|
| Cisplatin | SPL | Baseline |
| Transplatin | SPL | Trans isomer |
| cis-PtCl₂(en) | SPL | Bidentate ligand |
| Ferrocene | LIN | Haptic directionality |
| fac-Ir(ppy)₃ | OCT | `@b` axial tag |
| mer-Ir(ppy)₃ | OCT | `@b` axial tag |
| PdCl₂-R-BINAP | SPL | Axial chirality fixture |
| PdCl₂-RR-BDPP | SPL | P-chiral backbone |
| PdCl₂-RR-BDNN | SPL | N-chiral backbone |

## Acknowledgements

- **SCINE Molassembler**: 3D structure generation via distance geometry.
- **xyz2mol**: Jensen Group's algorithm for TMC graph generation from 3D coordinates.
- **RDKit**: Stereochemistry assignment and CIP oracle.
- **OpenBabel**: Chemical file handling.

## License

MIT
