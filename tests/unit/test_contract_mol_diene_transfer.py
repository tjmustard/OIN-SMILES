"""Regression: eta-alkene/diene/allyl bond-order transfer in ``build_contract_mol``.

``build_contract_mol`` recovers a generated ligand's bond orders by substructure-matching
``_flatten_template(template)`` into the generated fragment. BOTH sides are heavy-atom,
all-single connectivity graphs with equal atom counts, so the match is really an
*automorphism* search -- and RDKit returns an arbitrary one of them. For a symmetric eta
ligand most automorphisms move the bond orders: 1,5-cyclooctadiene's flattened 8-ring has
16 automorphisms and only 4 leave the two C=C on the metal-bound carbons. The other 12
rotate them onto the CH2-CH2 backbone, which is exactly what generated COD complexes
(GASBIN, PENGAT) re-encoded as: ``[CH]{1}1[CH]{1>}[CH2]=[CH2]...``. The eta3-bis-allyl
ABIRIO went further and put a double bond on a *methyl* (``C{3<}(=[CH3])``).

Note this is NOT the eta3-allyl failure of ``test_contract_mol_allyl_transfer`` -- there the
match returned EMPTY and nothing transferred. Here the match succeeds and the bond orders
land on the wrong atoms.

Two constraints pin the map, and both are needed:

* ``_oinSlot`` colouring -- a template atom binding OIN slot *s* may only map onto a
  generated atom binding that same slot. This kills the COD ring rotations, the para
  ring-swap that moved PIJCAO's alkyne onto a para-ethyl, and the macrocycle rotations of
  a porphyrin's four distinct N slots (a bare donor/non-donor colour leaves 8 of those
  legal and re-picking one shifts a slot label in the re-encoded OIN -- BOQPIG).
* a geometry score -- a *symmetric* eta3-allyl's three carbons all share one slot, so both
  slot-valid maps put the C=C on opposite ends of the allyl. The survivors are ranked
  against the embedded 3D geometry, where the true C=C is ~0.1 A shorter.

The slot colour additionally gates which template pairs with which fragment: equal
heavy-atom counts alone let a sigma-ethyl fragment consume an eta2-ethylene template.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from oinsmiles.generation.metallogen_adapter import (
    _flatten_template,
    _oin_fragment_templates,
    _select_match,
    _template_donor_slots,
    build_contract_mol,
)
from oinsmiles.generation.oin_parser import OINParser

# GASBIN-style eta2,eta2-COD on a square-planar Rh, with the NHC/Cl reduced to chlorides.
_COD_OIN = "[Rh_SPL].[Cl]{0}.[CH]{1}1=[CH]{1>}CC[CH]{2>}=[CH]{2}CC1.[Cl]{3}"
# ABAZEK-style eta3-allyl stripped of its phenyl, so the allyl is SYMMETRIC and no
# connectivity constraint can decide which end carries the C=C.
_BARE_ALLYL_OIN = "[Pd_TPL].[CH2]{0>}[CH]{0}=[CH]{0}.[Cl]{1}.[Cl]{2}"
# PIJCAO-style monodentate alkyne: the flat graph's para ring-swap maps the metal-bound
# alkyne carbon onto the terminal carbon of the para-ethyl.
_ALKYNE_OIN = "[Pt_SPL].C{0}#Cc1ccc(CC)cc1.[Cl]{1}.[Cl]{2}.[Cl]{3}"
# eta2-ethylene and a sigma-ethyl: same elements, same connectivity, same heavy-atom
# count -- only the donor count tells the two templates apart.
_ETHYLENE_ETHYL_OIN = "[Pt_SPL].[CH2]{0}=[CH2]{0}.[CH2]{1}C.[Cl]{2}.[Cl]{3}"


class _StubAtom:
    """Minimal stand-in for a MetalloGen atom (element + coordinate accessors)."""

    def __init__(self, element, coordinate):
        self._element = element
        self._coordinate = coordinate

    def get_element(self):
        return self._element

    def get_coordinate(self):
        return self._coordinate


class _StubMgMol:
    """Minimal stand-in for a MetalloGen result: atom_list + adj_matrix."""

    def __init__(self, atoms, bonds):
        self.atom_list = atoms
        n = len(atoms)
        adj = np.zeros((n, n), dtype=int)
        for i, j in bonds:
            adj[i, j] = adj[j, i] = 1
        self.adj_matrix = adj


def _template_for(parsed, n_atoms):
    """The parsed OIN's ligand template with ``n_atoms`` heavy atoms."""
    return next(t for t in _oin_fragment_templates(parsed) if t.GetNumAtoms() == n_atoms)


def _double_bond_pairs(mol):
    return {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in mol.GetBonds()
        if b.GetBondType() == Chem.BondType.DOUBLE
    }


def _cod_conformer():
    """A real 1,5-COD geometry; heavy-atom order is the ring order, C=C at (0,1) and (4,5)."""
    m = Chem.AddHs(Chem.MolFromSmiles("C1=CCCC=CCC1"))
    AllChem.EmbedMolecule(m, randomSeed=0xC0D)
    AllChem.MMFFOptimizeMolecule(m)
    m = Chem.RemoveHs(m)
    return m, m.GetConformer()


class TestCODAutomorphismPrecondition(unittest.TestCase):
    """Documents WHY the fix is needed: the flat COD query is 16-fold ambiguous."""

    def test_flat_cod_ring_has_bond_order_moving_automorphisms(self):
        parsed = OINParser().parse(_COD_OIN)
        t = _template_for(parsed, 8)
        doubles = {
            frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            for b in t.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE
        }
        self.assertEqual(len(doubles), 2, "eta2,eta2-COD template must carry two C=C")

        flat = _flatten_template(t)
        autos = flat.GetSubstructMatches(flat, uniquify=False, maxMatches=1000)
        self.assertEqual(len(autos), 16, "cyclooctane connectivity has a 16-element Aut group")

        keep = [a for a in autos if {frozenset(a[i] for i in e) for e in doubles} == doubles]
        self.assertEqual(
            len(keep), 4, "only 4 of the 16 automorphisms leave the C=C on the bound carbons"
        )
        # ...so an unconstrained GetSubstructMatch has a 12/16 chance of corrupting them.

    def test_donor_colour_alone_cannot_orient_a_symmetric_allyl(self):
        """Why the geometry score exists: both donor-valid maps of a bare allyl are legal."""
        parsed = OINParser().parse(_BARE_ALLYL_OIN)
        t = _template_for(parsed, 3)
        slots, _ = _template_donor_slots(t, parsed)
        self.assertEqual(set(slots), {0, 1, 2}, "all three allyl carbons bind the metal")
        self.assertEqual(set(slots.values()), {0}, "and they share one haptic slot")

        flat = _flatten_template(t, slots)
        qmol = Chem.MolFromSmiles("CCC")
        for a in qmol.GetAtoms():
            a.SetIntProp("_oinSlot", 0)
        params = Chem.SubstructMatchParameters()
        params.uniquify = False
        params.atomProperties = ["_oinSlot"]
        self.assertEqual(
            len(qmol.GetSubstructMatches(flat, params)),
            2,
            "the allyl reversal is slot-valid, so colouring leaves it ambiguous",
        )


class TestGeometrySelectsTheAutomorphism(unittest.TestCase):
    """``_select_match`` must rank donor-valid maps by the embedded 3D geometry."""

    def _allyl_pieces(self):
        parsed = OINParser().parse(_BARE_ALLYL_OIN)
        t = _template_for(parsed, 3)
        slots, _ = _template_donor_slots(t, parsed)
        flat = _flatten_template(t, slots)
        qmol = Chem.MolFromSmiles("CCC")
        for a in qmol.GetAtoms():
            a.SetIntProp("_oinSlot", 0)
        # template bonds: 0-1 single, 1=2 double
        return t, flat, qmol

    @staticmethod
    def _dmat(carr):
        return np.linalg.norm(carr[:, None, :] - carr[None, :, :], axis=-1)

    def test_double_bond_follows_the_short_bond(self):
        t, flat, qmol = self._allyl_pieces()
        q2g = [0, 1, 2]
        # C1=C2 short (1.34), C0-C1 long (1.50) -> the identity map is correct.
        carr = np.array([[0.0, 0.0, 0.0], [1.50, 0.0, 0.0], [2.84, 0.0, 0.0]])
        self.assertEqual(_select_match(qmol, flat, t, q2g, self._dmat(carr)), (0, 1, 2))

    def test_double_bond_follows_the_short_bond_when_reversed(self):
        t, flat, qmol = self._allyl_pieces()
        q2g = [0, 1, 2]
        # Mirror the geometry: now C0-C1 is the short bond, so the template's double
        # bond (1=2) must map onto generated (1,0) -- i.e. the reversed automorphism.
        carr = np.array([[0.0, 0.0, 0.0], [1.34, 0.0, 0.0], [2.84, 0.0, 0.0]])
        self.assertEqual(_select_match(qmol, flat, t, q2g, self._dmat(carr)), (2, 1, 0))

    def test_legacy_match_is_kept_when_only_noise_separates_the_candidates(self):
        """Repair, don't re-pick.

        Formal charges and _CIPCode stereo ride along on the chosen map, so re-picking
        among equally-good maps is observable downstream: a porphyrin's equivalent pyrrole
        rings admit many donor-valid maps whose scores differ only by FF noise, and
        swapping one shifted a donor's slot label in the re-encoded OIN (BOQPIG).
        """
        t, flat, qmol = self._allyl_pieces()
        q2g = [0, 1, 2]
        legacy = qmol.GetSubstructMatch(flat)
        # Both bonds ~equal: the reversed map is better by 0.01 A, far under SCORE_TOL.
        carr = np.array([[0.0, 0.0, 0.0], [1.40, 0.0, 0.0], [2.81, 0.0, 0.0]])
        self.assertEqual(
            _select_match(qmol, flat, t, q2g, self._dmat(carr)),
            legacy,
            "a sub-tolerance score gain must not churn the map the old code picked",
        )


class TestBuildContractMolTransfer(unittest.TestCase):
    """End-to-end through the real ``build_contract_mol`` with a stub MetalloGen result."""

    def test_cod_double_bonds_land_on_the_metal_bound_carbons(self):
        parsed = OINParser().parse(_COD_OIN)
        cod, conf = _cod_conformer()

        # Scramble the generated atom order relative to the template: generated carbon
        # (2 + perm[i]) corresponds to template ring atom i. Rotating by two is exactly
        # the corruption seen in GASBIN, so the identity map is NOT the answer.
        perm = [2, 3, 4, 5, 6, 7, 0, 1]
        pos = [np.array(conf.GetAtomPosition(i)) for i in range(8)]
        coords = [None] * 11
        for i in range(8):
            coords[2 + perm[i]] = pos[i]

        bound = [2 + perm[i] for i in (0, 1, 4, 5)]  # template donors -> generated indices
        centroid = np.mean([coords[i] for i in bound], axis=0)
        coords[0] = centroid + np.array([0.0, 0.0, 2.2])  # Rh
        coords[1] = coords[0] + np.array([2.4, 0.0, 0.0])  # Cl
        coords[10] = coords[0] + np.array([-2.4, 0.0, 0.0])  # Cl

        atoms = [_StubAtom("Rh", coords[0]), _StubAtom("Cl", coords[1])]
        atoms += [_StubAtom("C", coords[2 + k]) for k in range(8)]
        atoms.append(_StubAtom("Cl", coords[10]))

        ring = [(i, (i + 1) % 8) for i in range(8)]
        bonds = [(2 + perm[i], 2 + perm[j]) for i, j in ring]
        bonds += [(0, 1), (0, 10)] + [(0, b) for b in bound]

        mol = build_contract_mol(parsed, _StubMgMol(atoms, bonds))
        self.assertIsNotNone(mol, "contract mol must build")

        expected = {frozenset((2 + perm[0], 2 + perm[1])), frozenset((2 + perm[4], 2 + perm[5]))}
        self.assertEqual(
            _double_bond_pairs(mol),
            expected,
            "the two C=C must sit on the metal-bound carbon pairs, not the CH2-CH2 backbone",
        )
        for b in bound:
            self.assertEqual(
                mol.GetBondBetweenAtoms(0, b).GetBondType(),
                Chem.BondType.DATIVE,
                "metal->donor bonds stay dative",
            )

    def test_donor_count_separates_eta2_ethylene_from_sigma_ethyl(self):
        """Equal heavy-atom count + identical connectivity: only the donor count pairs them."""
        parsed = OINParser().parse(_ETHYLENE_ETHYL_OIN)
        # 0=Pt, 1,2=ethylene (both bound), 3,4=ethyl (only 3 bound), 5,6=Cl
        coords = [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 1.30, 2.05]),  # ethylene C (1.33 A apart)
            np.array([0.0, -0.03, 2.05]),
            np.array([2.05, 0.0, 0.0]),  # ethyl CH2
            np.array([3.58, 0.0, 0.0]),  # ethyl CH3 (1.53 A)
            np.array([-2.3, 0.0, 0.0]),
            np.array([0.0, 0.0, -2.3]),
        ]
        syms = ["Pt", "C", "C", "C", "C", "Cl", "Cl"]
        atoms = [_StubAtom(s, c) for s, c in zip(syms, coords)]
        bonds = [(1, 2), (3, 4), (0, 1), (0, 2), (0, 3), (0, 5), (0, 6)]

        mol = build_contract_mol(parsed, _StubMgMol(atoms, bonds))
        self.assertIsNotNone(mol)
        self.assertEqual(
            _double_bond_pairs(mol),
            {frozenset((1, 2))},
            "the C=C belongs to the two-donor ethylene, not the one-donor ethyl",
        )
        self.assertEqual(
            mol.GetBondBetweenAtoms(3, 4).GetBondType(),
            Chem.BondType.SINGLE,
            "the sigma-ethyl must not absorb the ethylene template's double bond",
        )


class TestMonodentateParaSwap(unittest.TestCase):
    """PIJCAO: a para ring-swap moved the alkyne onto the para-ethyl. Donor colour kills it."""

    def test_donor_colour_pins_the_alkyne_carbon(self):
        parsed = OINParser().parse(_ALKYNE_OIN)
        t = _template_for(parsed, 10)
        slots, _ = _template_donor_slots(t, parsed)
        self.assertEqual(set(slots), {0}, "only the terminal alkyne carbon binds the metal")

        flat = _flatten_template(t, slots)
        # The generated fragment: same skeleton, all single, slot colour on atom 0 only.
        qmol = _flatten_template(t, slots)
        uncoloured = qmol.GetSubstructMatches(flat, uniquify=False, maxMatches=100)
        params = Chem.SubstructMatchParameters()
        params.uniquify = False
        params.atomProperties = ["_oinSlot"]
        coloured = qmol.GetSubstructMatches(flat, params)

        self.assertGreater(
            len(uncoloured), len(coloured), "the para swap is a real automorphism of the flat graph"
        )
        self.assertTrue(coloured, "the identity map must survive the slot colour")
        for m in coloured:
            self.assertEqual(m[0], 0, "the alkyne carbon may only map onto a metal-binding atom")


if __name__ == "__main__":
    unittest.main()
