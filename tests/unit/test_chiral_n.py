"""N-chiral fixture tests: stability + RDKit CIP oracle.

Fixture: tests/fixtures/PdCl2-RR-BDNN.xyz  (64 atoms, K:1,2;)
  Complex: cis-[PdCl2(BDNN)] where BDNN = (2R,4R)-pentane-2,4-diyl-bis(diphenylamine)
  Chiral centres: C2 and C4 of the pentanediyl backbone (both R in IUPAC, S by CIP
  priority because of Ph substituents). N atoms are NOT CIP stereocentres — both carry
  two identical phenyl groups (and are tertiary amines without inversion barrier).

Candidate OIN (2026-03-04):
  [Pd_SPL].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}
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

# Human-reviewed expected OIN string (2026-03-04, generated from above XYZ)
_EXPECTED_OIN = (
    "[Pd_SPL].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}"
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
    def test_n_cip_oracle(self):
        """RDKit CIP oracle: backbone C stereocentres are both S.

        The N atoms in BDNN are tertiary amines with two identical phenyl
        groups and are not CIP stereocentres. The @/@@  encoding on the two
        backbone carbon atoms (C2, C4 of the pentanediyl chain) captures the
        chirality; both are assigned S by RDKit CIP.
        """
        oin = self._oin()
        ligand_smiles = extract_ligand_smiles(oin)

        mol = Chem.MolFromSmiles(ligand_smiles)
        self.assertIsNotNone(mol, f"RDKit could not parse: {ligand_smiles}")

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        chiral_codes = [
            a.GetPropsAsDict().get("_CIPCode")
            for a in mol.GetAtoms()
            if a.GetPropsAsDict().get("_CIPCode")
        ]
        # Both backbone carbons must be S (verified by RDKit CIP on 2026-03-04)
        self.assertEqual(
            sorted(chiral_codes),
            ["S", "S"],
            f"Expected [S, S] CIP codes, got {sorted(chiral_codes)} "
            f"from ligand SMILES: {ligand_smiles}",
        )


if __name__ == "__main__":
    unittest.main()
