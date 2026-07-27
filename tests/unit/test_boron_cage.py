"""The ``OIN_BORON_CAGE`` lever: deltahedral borane/carborane cages.

Background (measured, see ``docs/agentic-notes/v0.4.5/BORON_CAGE_v0.4.5.md``): the 34 `boron_cluster`
molecules in the v0.4.5 `encode_fail` cohort were documented as a permanent
ceiling on the grounds that RDKit has no 2c-2e Lewis structure for a 3c-2e cage.
That reasoning never applied. ``xyz2AC_obabel``'s covalent-radius distance
criterion perceives the cages *correctly*; a pruning loop then deletes an atom's
longest bonds while its connectivity exceeds ``max(atomic_valence[Z])``, and for
boron that cap is 4 against a cage vertex's 5-6 -- so 7-19 cage edges per molecule
are amputated and the cage shatters before any bond-order reasoning happens.

These tests pin both halves of the lever: OFF is inert, ON keeps the cage intact
and round-trips it.
"""

from __future__ import annotations

import os
import unittest
from pathlib import Path

from rdkit import Chem, RDLogger

from oinsmiles import XYZToSMILES
from oinsmiles.core.constants import TRANSITION_METALS_NUM
from oinsmiles.oin.compare import canonical_roundtrip_key
from oinsmiles.utils.aromaticity import OINEncodeError
from oinsmiles.utils.perception_core import boron_cage_vertices, read_xyz_file, xyz2AC_obabel

RDLogger.DisableLog("rdApp.*")

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"
#: closo-B12H12 dodecaborate on Rh -- all 30 cage edges are B-B.
OZAREO = FIXTURES / "OZAREO_comp_0.xyz"
#: o-carborane C2B10 on Ag(PPh3)2 -- a cage plus two ordinary phosphine ligands,
#: so it also pins that the relaxation does not leak into the non-cage fragments.
MODZUA = FIXTURES / "MODZUA_comp_0.xyz"
#: A Rh thiaborane (S is a cage vertex). Encoding it used to kill the process:
#: RDKit stamped a chiral tag on a 6-connected cage boron and `AssignStereochemistry`
#: then raised `RuntimeError: basic_string::_M_create` / aborted with
#: `free(): invalid size`. Pins `clear_boron_cage_stereo`.
KIXXOF = FIXTURES / "KIXXOF_comp_0.xyz"
#: A nido-C2B7 on Ru that is scored as a PASS today while being silently wrong:
#: with the lever off, 6 of its 12 B-B cage bonds are deleted and the encoder
#: invents a C=B double bond to balance the valences.
VEJXOZ = FIXTURES / "VEJXOZ_comp_0.xyz"

#: A boron-rich molecule with NO cage: an Ir boroxine, three borons in a B-O-B-O-B-O
#: ring, zero B-B bonds. The end-to-end control for the motif gate.
ASUVIV = FIXTURES / "ASUVIV_comp_0.xyz"
#: An Fe complex carrying [BH3-] / [BH-] / borane groups -- boron at several
#: coordination numbers, still no deltahedron.
AROTAE = FIXTURES / "AROTAE_comp_0.xyz"

LEVER = "OIN_BORON_CAGE"


class _LeverMixin(unittest.TestCase):
    """Sets ``OIN_BORON_CAGE`` EXPLICITLY in both directions.

    Deleting the variable used to mean "off" and stopped meaning that when the lever was
    promoted to default-ON in v0.4.6: every ``test_lever_off_*`` below silently became a second
    lever-ON test and asserted the amputated-cage behaviour against the fixed path. Third and
    fourth occurrences of that trap in this release -- see
    ``test_levers.py::TestNoTestUnsetsAPromotedLever``, which now makes it a lint.
    """

    def setUp(self):
        self._saved = os.environ.get(LEVER)
        os.environ[LEVER] = "0"

    def tearDown(self):
        if self._saved is None:
            os.environ.pop(LEVER, None)
        else:
            os.environ[LEVER] = self._saved

    def set_lever(self, on):
        os.environ[LEVER] = "1" if on else "0"


class TestCageMotifDetection(unittest.TestCase):
    """``boron_cage_vertices`` recognises the motif, not merely "has boron"."""

    def _vertices(self, smiles):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
        return boron_cage_vertices(atoms, Chem.rdmolops.GetAdjacencyMatrix(mol))

    def test_borate_anion_is_not_a_cage(self):
        # BPh4-: four B-C bonds, no B-B bond at all.
        self.assertEqual(self._vertices("[B-](c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1"), set())

    def test_tetrafluoroborate_is_not_a_cage(self):
        self.assertEqual(self._vertices("[B-](F)(F)(F)F"), set())

    def test_diboron_bond_is_not_a_cage(self):
        # A B-B bond with no third boron cannot close a triangle.
        self.assertEqual(self._vertices("OB(O)B(O)O"), set())

    def test_boroxine_ring_is_not_a_cage(self):
        # B-O-B alternation: three borons, zero B-B bonds.
        self.assertEqual(self._vertices("B1OB(O)OB(O)O1"), set())

    def test_linear_three_boron_chain_is_not_a_cage(self):
        # Three borons, two B-B bonds, but no triangle -- stricter than the
        # pre-existing `_is_electron_deficient_cluster` (>=3 B and >=1 B-B).
        self.assertEqual(self._vertices("OB(O)B(O)B(O)O"), set())

    def test_boron_triangle_is_a_cage(self):
        verts = self._vertices("B1B(O)B1O")
        self.assertEqual(len(verts), 3)


class TestPruningExemption(_LeverMixin):
    """The pruning loop is what destroys the cage; the lever is what stops it."""

    def _b_stats(self, path):
        """Max boron degree and B-B edge count in the encoder's own adjacency matrix.

        The metal is excluded from the degree so the number is a property of the
        cage, not of how many vertices happen to face the metal.
        """
        atoms, _charge, xyz = read_xyz_file(str(path))
        AC, _mol = xyz2AC_obabel(atoms, xyz, tolerance=0.5)
        metals = {i for i, z in enumerate(atoms) if z in TRANSITION_METALS_NUM}
        b = [i for i, z in enumerate(atoms) if z == 5]
        max_deg = max(sum(int(AC[i][j]) for j in range(len(atoms)) if j not in metals) for i in b)
        n_bb = sum(int(AC[i][j]) for k, i in enumerate(b) for j in b[k + 1 :])
        return max_deg, n_bb

    def test_lever_off_amputates_the_cage(self):
        self.set_lever(False)
        max_deg, n_bb = self._b_stats(OZAREO)
        # Every boron forced down to the atomic_valence[5] cap of 4...
        self.assertEqual(max_deg, 4)
        # ...costing more than a third of the 30 icosahedral B-B edges.
        self.assertLess(n_bb, 30)

    def test_lever_on_keeps_the_full_icosahedron(self):
        self.set_lever(True)
        max_deg, n_bb = self._b_stats(OZAREO)
        # closo-B12H12: each vertex has 5 cage neighbours plus one exo substituent
        # (an H, or the amide N on the single substituted vertex).
        self.assertEqual(max_deg, 6)
        # A B12 icosahedron has exactly 30 edges, all of them B-B.
        self.assertEqual(n_bb, 30)


class TestCageEncodes(_LeverMixin):
    def test_lever_off_still_raises_the_typed_error(self):
        self.set_lever(False)
        with self.assertRaises(OINEncodeError):
            XYZToSMILES().convert(str(OZAREO))

    def test_lever_on_encodes_and_is_deterministic(self):
        self.set_lever(True)
        first = XYZToSMILES().convert(str(OZAREO))
        self.assertTrue(first)
        self.assertEqual(first, XYZToSMILES().convert(str(OZAREO)))

    def test_encoded_cage_keeps_all_twelve_borons(self):
        self.set_lever(True)
        oin = XYZToSMILES().convert(str(OZAREO))
        # 12 cage borons, written as `B`/`[BH]`; count bracketed and bare alike.
        n_b = sum(1 for i, ch in enumerate(oin) if ch == "B" and not oin[i:].startswith("Br"))
        self.assertEqual(n_b, 12)

    def test_key_does_not_degrade_to_the_raw_fallback(self):
        # A cage fragment that cannot be parsed contributes its literal input
        # SMILES to the key, which makes the key atom-order dependent. The lever's
        # `_parse_fragment` rung exists to prevent exactly that.
        self.set_lever(True)
        oin = XYZToSMILES().convert(str(OZAREO))
        key = canonical_roundtrip_key(oin)
        self.assertNotIn("RAW:", str(key))
        self.assertEqual(str(key), str(canonical_roundtrip_key(oin)))

    def test_mixed_complex_encodes_cage_and_ordinary_ligands(self):
        # MODZUA is Ag(PPh3)2 + a carborane: the two phosphines must still come
        # out as ordinary aromatic fragments, so the relaxation is not leaking.
        self.set_lever(True)
        oin = XYZToSMILES().convert(str(MODZUA))
        self.assertIn("c1ccc(P", oin)
        self.assertIn("[BH]", oin)


class TestCageStereoMustNeverBeTagged(_LeverMixin):
    """A chiral tag on a cage vertex is a native memory fault, not a bad descriptor.

    Measured before ``clear_boron_cage_stereo`` existed: ``KIXXOF`` raised
    ``RuntimeError: basic_string::_M_create`` out of ``Chem.AssignStereochemistry``,
    and its sibling ``DUDTIG`` killed the interpreter with ``free(): invalid size``
    (SIGABRT). RDKit has no stereo permutation table for a 5-/6-connected boron.
    """

    def test_thiaborane_encodes_without_crashing(self):
        self.set_lever(True)
        oin = XYZToSMILES().convert(str(KIXXOF))
        self.assertTrue(oin)

    def test_no_chiral_tag_survives_on_a_cage_vertex(self):
        self.set_lever(True)
        oin = XYZToSMILES().convert(str(KIXXOF))
        # `[B@]` / `[B@@]` / `[B@H]` would mean a tag reached the serializer.
        self.assertNotIn("[B@", oin)


class TestSilentCorruptionOfAPassingMolecule(_LeverMixin):
    """The pruning defect also corrupts molecules that are scored as passing.

    ``VEJXOZ`` produces an OIN either way and round-trips against its own mol
    either way -- which is exactly why the round-trip key could not see that the
    lever-off graph is wrong. The check has to be on the cage bond count.
    """

    def _cage_bond_count(self, path):
        from oinsmiles.utils.perception_tmc import get_tmc_mol

        mol, _xyz = get_tmc_mol(str(path), 0, with_stereo=False)
        bb = sum(
            1
            for b in mol.GetBonds()
            if b.GetBeginAtom().GetAtomicNum() == 5 and b.GetEndAtom().GetAtomicNum() == 5
        )
        spurious = sum(
            1
            for b in mol.GetBonds()
            if str(b.GetBondType()) == "DOUBLE"
            and 5 in (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum())
        )
        return bb, spurious

    def test_lever_off_deletes_half_the_cage_and_invents_a_double_bond(self):
        self.set_lever(False)
        bb, spurious = self._cage_bond_count(VEJXOZ)
        self.assertEqual(bb, 6)  # geometry has 12
        self.assertGreaterEqual(spurious, 1)  # a C=B double bond in a carborane

    def test_lever_on_keeps_the_cage_and_invents_nothing(self):
        self.set_lever(True)
        bb, spurious = self._cage_bond_count(VEJXOZ)
        self.assertEqual(bb, 12)
        self.assertEqual(spurious, 0)


class TestNonCageMoleculesUnaffected(_LeverMixin):
    """The lever must be inert for anything without the motif, ON or OFF.

    Four golden fixtures, each re-encoded in both arms and compared to itself --
    the corpus-scale version of this check is ``tools/boron_regression_ab.py``.
    """

    GOLDENS = ["CisPlatin.xyz", "Ferrocene.xyz", "fac-Ir(ppy)3.xyz", "Cis-PtCl2(en).xyz"]

    def _both_arms(self, path):
        self.set_lever(False)
        off = XYZToSMILES().convert(str(path))
        self.set_lever(True)
        on = XYZToSMILES().convert(str(path))
        return off, on

    def test_goldens_byte_identical_across_the_lever(self):
        for name in self.GOLDENS:
            off, on = self._both_arms(FIXTURES / name)
            self.assertEqual(off, on, f"{name} changed when OIN_BORON_CAGE was set")

    def test_boron_rich_non_cage_molecules_byte_identical(self):
        """The sharpest control: boron everywhere, no deltahedron, must not move.

        A golden fixture has no boron at all, so it cannot show whether the gate is
        the *motif* or merely the *element*. These two can. ASUVIV is a boroxine
        (three borons, zero B-B bonds); AROTAE carries borohydride/borane groups at
        several coordination numbers. Both are byte-identical across the lever, which
        is what makes "scoped to the B-B-B triangle" a measurement rather than a
        claim about the code.
        """
        for path in (ASUVIV, AROTAE):
            off, on = self._both_arms(path)
            self.assertTrue(off, f"{path.name} did not encode at all")
            self.assertEqual(off, on, f"{path.name} changed when OIN_BORON_CAGE was set")


if __name__ == "__main__":
    unittest.main()


class TestValenceBypassIsScopedToBoron(_LeverMixin):
    """The bypass must not leak to non-boron fragments. Regression guard for a real leak.

    ``_parse_fragment``'s cage rung skips ``SANITIZE_PROPERTIES`` -- a valence-RULE bypass. When
    ``OIN_BORON_CAGE`` was promoted to default-ON it was unscoped, so it applied to every fragment
    and parsed ``C#O``: carbon monoxide fails the valence check and nothing else. CO is among the
    commonest ligands in transition-metal chemistry, so the promotion was silently changing
    chemistry far outside the 34 boron molecules it was justified on.
    """

    def test_carbon_monoxide_still_falls_back_with_the_lever_ON(self):
        from oinsmiles.oin.compare import _parse_fragment

        self.set_lever(True)
        self.assertIsNone(
            _parse_fragment("C#O"),
            "over-valent C#O must NOT take the boron valence bypass -- it contains no boron, so "
            "it cannot be a deltahedral cage, and parsing it hides an over-valent carbon that "
            "the RAW: fallback exists to exclude",
        )

    # The other half -- that scoping costs the boron population nothing -- is already covered by
    # TestPruningExemption::test_lever_on_keeps_the_full_icosahedron (30/30 B-B edges survive) and
    # TestCageEncodes::test_key_does_not_degrade_to_the_raw_fallback, both of which exercise a
    # REAL cage fixture with the lever on. A hand-written cage SMILES was tried here and only
    # proved that hand-written cage SMILES are easy to get wrong.
