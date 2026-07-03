"""Axial-chirality handling tests for BINAP-containing Pd complex.

Per MiniPRD_ChiralTests (US-004 / Test 5): axial chirality descriptors are a
Novel output requiring human review before hard assertions can be written.
The candidate artifact is saved to tests/candidate_outputs/axial_chiral_encoded.smi.

The current OIN pipeline does NOT encode biaryl atropisomeric chirality (no
`/` or `\\` bond-stereo descriptors in the BINAP substructure).  The skip
decorators will be removed and replaced with concrete assertions once a
human has reviewed the candidate artifact and the pipeline encodes axial
chirality.

Fixture: tests/fixtures/PdCl2-R-BINAP.xyz  (81 atoms, K:19,38;, R-BINAP)
Candidate artifact: tests/candidate_outputs/axial_chiral_encoded.smi (2026-03-04)
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from oinsmiles import XYZToSMILES
from tests.unit.test_helpers import extract_ligand_smiles

try:
    from rdkit import Chem
    from rdkit.Chem import FindPotentialStereo, StereoInfo

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

_BINAP_XYZ = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../fixtures/PdCl2-R-BINAP.xyz")
)

_CANDIDATE_ARTIFACT = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__), "../../tests/candidate_outputs/axial_chiral_encoded.smi"
    )
)


class TestAxialChiral(unittest.TestCase):
    def _oin(self) -> str:
        return XYZToSMILES().convert(_BINAP_XYZ)

    def test_pipeline_produces_output(self):
        """Pipeline must return a non-None OIN string (non-crash guard)."""
        result = self._oin()
        self.assertIsNotNone(result)
        self.assertIsInstance(result, str)

    def test_candidate_artifact_saved(self):
        """Candidate artifact file must exist for human review."""
        self.assertTrue(
            os.path.exists(_CANDIDATE_ARTIFACT),
            f"Candidate artifact not found: {_CANDIDATE_ARTIFACT}",
        )

    @unittest.skip(
        "Axial chirality not yet encoded by pipeline. "
        "Un-skip after pipeline update and human review of "
        "tests/candidate_outputs/axial_chiral_encoded.smi."
    )
    @unittest.skipUnless(_RDKIT_AVAILABLE, "rdkit not installed")
    def test_axial_chiral_descriptor_present(self):
        """OIN ligand SMILES must contain axial-chirality bond descriptors.

        BINAP has an R-configured biaryl atropisomeric axis (K:19,38 in the
        XYZ comment).  Once the pipeline encodes this, the ligand SMILES must
        contain `/` or `\\` bond stereo tokens at the biaryl bond, and
        Chem.FindPotentialStereo() must identify at least one atropisomeric
        stereo element.
        """
        oin = self._oin()
        ligand_smiles = extract_ligand_smiles(oin)

        # Axial chirality requires bond-stereo tokens
        self.assertTrue(
            "/" in ligand_smiles or "\\" in ligand_smiles,
            f"No bond-stereo token in ligand SMILES: {ligand_smiles}",
        )

        mol = Chem.MolFromSmiles(ligand_smiles)
        self.assertIsNotNone(mol)
        stereo_info = FindPotentialStereo(mol)
        atropisomeric = [
            s for s in stereo_info if s.type == StereoInfo.StereoType.Bond_Atropisomeric
        ]
        self.assertGreater(
            len(atropisomeric),
            0,
            "No atropisomeric stereo element found in BINAP ligand SMILES.",
        )


if __name__ == "__main__":
    unittest.main()
