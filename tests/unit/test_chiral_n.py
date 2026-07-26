"""N-chiral fixture tests: stability + a CIP check anchored on the 3D geometry.

Fixture: tests/fixtures/PdCl2-RR-BDNN.xyz  (64 atoms, K:1,2;)
  Complex: cis-[PdCl2(BDNN)] where BDNN = (2R,4R)-pentane-2,4-diyl-bis(diphenylamine)
  Chiral centres: C2 and C4 of the pentanediyl backbone, **both R**. N atoms are NOT CIP
  stereocentres — both carry two identical phenyl groups (and are tertiary amines without an
  inversion barrier).

⚠ Same correction as ``test_chiral_p`` — read that module's docstring for the full account.
In short: the old ``["S", "S"]`` assertion ran CIP on a SMILES reparsed from the encoder's own
output, which relabels a tag rather than checking it, so an inverted tag passed. The geometry
is the arbiter and says R at both centres, agreeing with the ``(2R,4R)`` in the filename.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES
from tests.unit.test_helpers import extract_ligand_smiles

try:
    from rdkit import Chem

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

_BDNN_XYZ = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../fixtures/PdCl2-RR-BDNN.xyz")
)

# Expected OIN. Updated 2026-07-26 (v0.4.5 Lane 8): both backbone tags flipped when
# OIN_STABLE_STEREO was promoted to default-ON; the previous golden encoded the inverted
# configuration. Both spellings describe R at both centres -- the per-atom symbol differs from
# the label because parity is relative to neighbour write order.
_EXPECTED_OIN = (
    "[Pd_SPL].C[C@H](C[C@@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}"
)


class TestChiralN(unittest.TestCase):
    def _oin(self) -> str:
        return XYZToSMILES().convert(_BDNN_XYZ)

    def test_n_stability(self):
        """Round-trip: encode XYZ → OIN matches human-reviewed candidate."""
        self.assertEqual(self._oin(), _EXPECTED_OIN)

    def test_n_no_spurious_chiral_tag(self):
        """Negative control (Stereo Phase 4): nitrogen is explicitly out of
        scope for the Zone-A lone-pair CIP feature (trivalent [N@] is
        RDKit-cleared as non-stereogenic amine inversion) -- Zone-A N must
        keep today's clearing behaviour. Byte-identical golden (above) plus
        this explicit tag-absence assertion.
        """
        oin = self._oin()
        self.assertNotIn("[N@]", oin)
        self.assertNotIn("[N@@]", oin)

    @unittest.skipUnless(_RDKIT_AVAILABLE, "rdkit not installed")
    def test_n_cip_from_geometry(self):
        """Ground truth, from coordinates: both backbone C stereocentres are R.

        The N atoms in BDNN are tertiary amines with two identical phenyl groups and are not
        CIP stereocentres; the @/@@ on the two backbone carbons carries the configuration.
        Labelled with ``rdCIPLabeler`` from ``AssignStereochemistryFrom3D``, independent of the
        emitted string -- see the module docstring.
        """
        from rdkit.Chem import rdCIPLabeler

        from oinsmiles.utils.xyz2mol import get_tmc_mol

        parent = get_tmc_mol(_BDNN_XYZ, 0, with_stereo=True)
        if isinstance(parent, tuple):
            parent = parent[0]
        Chem.AssignStereochemistryFrom3D(parent)
        rdCIPLabeler.AssignCIPLabels(parent)
        chiral_codes = [
            a.GetPropsAsDict()["_CIPCode"]
            for a in parent.GetAtoms()
            if a.HasProp("_CIPCode") and a.GetSymbol() == "C"
        ]
        self.assertEqual(
            sorted(chiral_codes),
            ["R", "R"],
            f"Expected [R, R] from the 3D geometry, got {sorted(chiral_codes)}. "
            "The fixture is (2R,4R); if this fails the fixture or the perception changed, "
            "NOT the golden string -- do not 'fix' it by editing the expectation.",
        )

    @unittest.skipUnless(_RDKIT_AVAILABLE, "rdkit not installed")
    def test_emitted_string_agrees_with_the_geometry(self):
        """The emitted tags must denote the configuration the coordinates do."""
        from rdkit.Chem import rdCIPLabeler

        mol = Chem.MolFromSmiles(extract_ligand_smiles(self._oin()))
        self.assertIsNotNone(mol, "emitted ligand body must be parseable")
        rdCIPLabeler.AssignCIPLabels(mol)
        emitted = sorted(
            a.GetPropsAsDict()["_CIPCode"]
            for a in mol.GetAtoms()
            if a.HasProp("_CIPCode") and a.GetSymbol() == "C"
        )
        self.assertEqual(
            emitted, ["R", "R"], f"emitted string says {emitted}, geometry says [R, R]"
        )


if __name__ == "__main__":
    unittest.main()
