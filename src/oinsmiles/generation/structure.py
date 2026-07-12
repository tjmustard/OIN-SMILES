"""Shared result type for the OIN 3D-generation backends.

Kept in a lightweight module (no heavy engine imports) so that ``engine.py``
can reference the return type at module load time while still importing the
concrete ``MetalloGenAdapter`` lazily inside ``OIN3DGenerator.__init__``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from rdkit import Chem


@dataclass
class GeneratedStructure:
    """Result of 3D structure generation: XYZ block + optional bonded RDKit mol."""

    xyz: str
    mol: Optional[Chem.Mol] = None
