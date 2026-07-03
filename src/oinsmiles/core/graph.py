from enum import Enum, auto
from typing import Dict, List, Optional, Tuple

import numpy as np


class BondType(Enum):
    """Enumeration of bond types in a transition metal complex graph."""

    COVALENT = auto()
    DATIVE = auto()
    ZERO_ORDER = auto()


class Atom:
    """A single atom in a TMC graph, with its coordinates and bonds."""

    def __init__(self, index: int, symbol: str, coords: Tuple[float, float, float]):
        """Initialize an atom with its index, element symbol, and 3D coordinates."""
        self.index = index
        self.symbol = symbol
        self.coords = coords
        self.neighbors: List["Atom"] = []
        self.bonds: Dict[int, BondType] = {}  # Map neighbor index to BondType

    def add_neighbor(self, neighbor: "Atom", bond_type: BondType):
        """Record a bond of the given type from this atom to a neighbor."""
        self.neighbors.append(neighbor)
        self.bonds[neighbor.index] = bond_type


class TMCGraph:
    """Molecular graph of a transition metal complex."""

    def __init__(self):
        """Initialize an empty graph with no atoms and no metal center."""
        self.atoms: List[Atom] = []
        self.metal_index: Optional[int] = None
        self.relative_coords: Dict[
            int, Tuple[float, float, float]
        ] = {}  # Map atom index to relative coords

    def add_atom(self, symbol: str, coords: Tuple[float, float, float]) -> int:
        """Add an atom and return its assigned index."""
        index = len(self.atoms)
        atom = Atom(index, symbol, coords)
        self.atoms.append(atom)
        return index

    def add_bond(self, idx1: int, idx2: int, bond_type: BondType = BondType.COVALENT):
        """Add a bidirectional bond of the given type between two atoms."""
        atom1 = self.atoms[idx1]
        atom2 = self.atoms[idx2]
        atom1.add_neighbor(atom2, bond_type)
        atom2.add_neighbor(atom1, bond_type)

    def set_metal_center(self, index: int):
        """Mark the atom at the given index as the metal center."""
        self.metal_index = index

    def calculate_relative_coords(self):
        """Calculates coordinates relative to the geometric center (centroid)."""
        if not self.atoms:
            return

        # Calculate centroid
        coords = np.array([atom.coords for atom in self.atoms])
        centroid = np.mean(coords, axis=0)

        for atom in self.atoms:
            rel_vec = np.array(atom.coords) - centroid
            self.relative_coords[atom.index] = tuple(rel_vec)

    def get_smiles(self) -> str:
        """Placeholder for SMILES generation.

        In a real implementation, this would use RDKit or OpenBabel.
        """
        # TODO: Integrate with RDKit to generate canonical SMILES from the graph
        return ""
