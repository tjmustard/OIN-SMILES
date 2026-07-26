"""P-chiral fixture tests: stability + a CIP check anchored on the 3D geometry.

Fixture: tests/fixtures/PdCl2-RR-BDPP.xyz  (64 atoms, K:1,2;)
  Complex: cis-[PdCl2(BDPP)] where BDPP = (2R,4R)-pentane-2,4-diyl-bis(diphenylphosphino)
  Chiral centres: C2 and C4 of the pentanediyl backbone, **both R**. P atoms are NOT CIP
  stereocentres — both carry two identical phenyl groups.

⚠ HISTORY — this file pinned the wrong answer for four months, and how it did so is worth
keeping. Until v0.4.5 it asserted ``["S", "S"]``, citing "verified by RDKit CIP". The
verification was circular: it took the encoder's OWN emitted string, reparsed it with
``MolFromSmiles``, and ran CIP on the result. ``rdCIPLabeler`` converts a parity tag into an
R/S label — it does not check that tag against anything. Hand it an inverted tag and it
returns an inverted label with full confidence. So the "oracle" was a snapshot of the
encoder's output, and when Lane 8 fixed the underlying tag instability the snapshot is what
failed.

The arbiter has to be the geometry, which is the one thing no encoder bug can rewrite:
``AssignStereochemistryFrom3D`` on the parent complex gives R at both centres, agreeing with
the ``(2R,4R)`` in the fixture's own name. ``test_p_cip_from_geometry`` below now does that,
and the string assertion is a consistency check against it rather than the primary claim.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES
from tests.unit.test_helpers import extract_ligand_smiles

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem  # noqa: F401  (import probe for RDKit availability)

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

_BDPP_XYZ = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../fixtures/PdCl2-RR-BDPP.xyz")
)

# Expected OIN. Updated 2026-07-26 (v0.4.5 Lane 8): both backbone tags flipped, @@/@ -> @/@@,
# when OIN_STABLE_STEREO was promoted to default-ON. The previous golden encoded the inverted
# configuration -- see the module docstring. Both spellings describe R at both centres; the
# per-atom symbol differs from the label because parity is relative to neighbour write order.
_EXPECTED_OIN = (
    "[Pd_SPL].C[C@H](C[C@@H](C)P{0}(c1ccccc1)c1ccccc1)P{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}"
)


class TestChiralP(unittest.TestCase):
    def _oin(self) -> str:
        return XYZToSMILES().convert(_BDPP_XYZ)

    def test_p_stability(self):
        """Round-trip: encode XYZ → OIN matches human-reviewed candidate."""
        self.assertEqual(self._oin(), _EXPECTED_OIN)

    def test_p_no_spurious_zone_a_tag(self):
        """Negative control (Stereo Phase 4, Q4/RISK-4): BDPP's P atoms bind
        Pd directly (Zone-A) but are NOT CIP stereocentres (two identical
        phenyl groups). The Zone-A P lone-pair CIP feature must not leak a
        spurious [P@]/[P@@] tag onto them -- byte-identical golden (above)
        plus this explicit tag-absence assertion, per the MiniPRD's negative
        constraints.
        """
        oin = self._oin()
        self.assertNotIn("[P@]", oin)
        self.assertNotIn("[P@@]", oin)

    @unittest.skipUnless(_RDKIT_AVAILABLE, "rdkit not installed")
    def test_p_cip_from_geometry(self):
        """Ground truth, from coordinates: both backbone C stereocentres are R.

        Derived with ``AssignStereochemistryFrom3D`` on the parent complex and labelled with
        ``rdCIPLabeler`` (the rigorous implementation, not the legacy ``AssignStereochemistry``).
        This is deliberately independent of anything the encoder emits -- see the module
        docstring for the four months this file spent asserting the encoder against itself.
        Agrees with the ``(2R,4R)`` in the fixture's filename.
        """
        from rdkit.Chem import rdCIPLabeler

        from oinsmiles.utils.xyz2mol import get_tmc_mol

        parent = get_tmc_mol(_BDPP_XYZ, 0, with_stereo=True)
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
        """The emitted tags must denote the same configuration the coordinates do.

        Closes the loop the old circular oracle left open: it checked the string against
        itself, so an inverted tag was self-consistent and passed. Here the string is read
        back and required to agree with the independent geometric answer above.
        """
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
