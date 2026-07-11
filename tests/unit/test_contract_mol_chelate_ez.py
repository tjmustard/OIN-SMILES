"""Regression: ``build_contract_mol`` must drop E/Z from metal-ring-locked double bonds.

The fast re-encode feeds ``build_contract_mol``'s output straight into
``get_oin_string``. ``AssignStereochemistryFrom3D`` (run inside build_contract_mol)
stamps a directional E/Z marker on every localized ring double bond -- including a
porphyrin's meso ``C=C``/``C=N`` bridges and any chelate/eta-locked alkene, which are
physically ring-locked and have no free E/Z. The XYZ->OIN *forward* encode strips those
(``translator.XYZToSMILES.convert`` -> ``_clear_chelate_locked_bond_stereo``), so if the
generator path does not, the round trip fails on a bond that was never stereogenic --
which is exactly what the ``macrocycle_perception`` E/Z-only rows (ABIPIK, ESANAO,
IFIYIH, JIVNOR, LOLROW, QAKLET, ...) re-encoded with.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np
from rdkit import Chem

from oinsmiles.generation.metallogen_adapter import build_contract_mol
from oinsmiles.generation.oin_parser import OINParser

# A 1,4-diaza-1,3-butadiene (DAD / di-imine) chelate on square-planar Ni, plus two
# chlorides. Both imine N bind the metal, so the backbone N=C-C=N closes a 5-membered
# ring THROUGH the metal and its two C=N are ring-locked.
_DAD_OIN = "[Ni_SPL].CN{0}=CC=N{1}C.[Cl]{2}.[Cl]{3}"


class _StubAtom:
    def __init__(self, element, coordinate):
        self._element = element
        self._coordinate = coordinate

    def get_element(self):
        return self._element

    def get_coordinate(self):
        return self._coordinate


class _StubMgMol:
    def __init__(self, atoms, bonds):
        self.atom_list = atoms
        n = len(atoms)
        adj = np.zeros((n, n), dtype=int)
        for i, j in bonds:
            adj[i, j] = adj[j, i] = 1
        self.adj_matrix = adj


def _dad_stub():
    """Ni-DAD-Cl2 with a deliberately puckered backbone so 3D perception WOULD assign
    E/Z to the two ring-locked C=N bonds (which the fix must then clear)."""
    # idx: 0=Ni 1=N 2=C 3=C 4=N 5=Cme(on N1) 6=Cme(on N4) 7=Cl 8=Cl
    coords = {
        0: (0.0, 0.0, 0.0),  # Ni
        1: (1.15, 1.15, 0.10),  # N
        2: (2.45, 0.75, -0.35),  # C (=N1)
        3: (2.45, -0.75, -0.35),  # C (=N4), backbone puckered in z
        4: (1.15, -1.15, 0.30),  # N
        5: (0.95, 2.55, 0.55),  # methyl on N1
        6: (0.95, -2.55, -0.55),  # methyl on N4
        7: (-1.15, 1.15, 0.0),  # Cl
        8: (-1.15, -1.15, 0.0),  # Cl
    }
    atoms = [
        _StubAtom("Ni", np.array(coords[0])),
        _StubAtom("N", np.array(coords[1])),
        _StubAtom("C", np.array(coords[2])),
        _StubAtom("C", np.array(coords[3])),
        _StubAtom("N", np.array(coords[4])),
        _StubAtom("C", np.array(coords[5])),
        _StubAtom("C", np.array(coords[6])),
        _StubAtom("Cl", np.array(coords[7])),
        _StubAtom("Cl", np.array(coords[8])),
    ]
    bonds = [
        (1, 2),
        (2, 3),
        (3, 4),  # N=C-C=N backbone
        (1, 5),
        (4, 6),  # N-methyls
        (0, 1),
        (0, 4),
        (0, 7),
        (0, 8),  # metal bonds
    ]
    return _StubMgMol(atoms, bonds)


def _ring_double_bond_stereo(mol):
    """(bond_pair -> stereo) for every C=N/C=C double bond in a metal ring."""
    out = {}
    for b in mol.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            out[frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))] = b.GetStereo()
    return out


class TestChelateLockedEZ(unittest.TestCase):
    def test_metal_ring_locked_double_bonds_carry_no_ez(self):
        parsed = OINParser().parse(_DAD_OIN)
        mol = build_contract_mol(parsed, _dad_stub())
        self.assertIsNotNone(mol, "contract mol must build")

        stereo = _ring_double_bond_stereo(mol)
        # The backbone must have been perceived as two C=N double bonds.
        self.assertEqual(len(stereo), 2, f"expected two C=N double bonds, got {stereo}")
        for pair, st in stereo.items():
            self.assertEqual(
                st,
                Chem.BondStereo.STEREONONE,
                f"ring-locked C=N {set(pair)} must have E/Z cleared, got {st}",
            )


if __name__ == "__main__":
    unittest.main()
