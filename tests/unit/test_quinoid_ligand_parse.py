"""Unit tests for MetalloGen parsing of un-kekulizable aromatic ligands.

Ti/Hf amidinate & 2-iminopyridine complexes encode a ligand whose aromatic
ring carries an exocyclic ``C=N`` (a quinoid system with NO valid aromatic
Kekule structure). ``process.get_ace_mol_from_rd_mol`` could not kekulize it,
reverted to the aromatic mol, then raised ``KeyError:
rdkit.Chem.rdchem.BondType.AROMATIC`` when reading bond orders -- surfacing as
"MetalloGen failed to generate any conformers" (the Track C ``no_conformers``
bucket: ABERIK/ABEROQ/ABERUW/ACUYUU/AFAMEB/AFAMIF/AFAMOL/AFIXUJ).

The fix approximates un-kekulizable aromatic bonds as single order instead of
crashing; the explicit exocyclic ``C=N`` is preserved and the final re-encoded
OIN is perceived from the generated 3D geometry, not this bond-order matrix.
"""

import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

from rdkit import Chem

from oinsmiles.generator3d import om, process

# The offending quinoid ligand (2-iminopyridine): aromatic ring + exocyclic C=N.
_QUINOID_LIG = "Cc1cc(C)c(-c2cccc(=[N:5]c3c(C)cc(C)cc3C)[n:3]2)c(C)c1"
# ABERIK's real m-SMILES, which raised the KeyError before the fix.
_ABERIK_MSMILES = (
    "[Ti]|[Cl:1]|[Cl:2]|[Cl:4]|"
    "Cc1cc(C)c(-c2cccc(=[N:5]c3c(C)cc(C)cc3C)[n:3]2)c(C)c1|"
    "C[NH:6]C|6_octahedral"
)


class TestGetAceMolAromaticFallback(unittest.TestCase):
    def test_unkekulizable_aromatic_does_not_raise(self):
        # Reproduce the exact input state get_ligand_from_smiles builds.
        mol = Chem.MolFromSmiles(_QUINOID_LIG, sanitize=False)
        self.assertIsNotNone(mol)
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.AddHs(mol, explicitOnly=False, addCoords=False)
        # Before the fix this raised KeyError(BondType.AROMATIC).
        ace = process.get_ace_mol_from_rd_mol(mol)
        self.assertGreater(len(ace.atom_list), 0)
        # A bond-order matrix was produced for every atom.
        self.assertEqual(ace.bo_matrix.shape[0], len(ace.atom_list))


class TestOmParseQuinoid(unittest.TestCase):
    def test_full_msmiles_parses_without_keyerror(self):
        # This is the call inside generate_3d_structures that raised, causing the
        # adapter to report "failed to generate any conformers".
        metal_complex = om.get_om_from_modified_smiles(_ABERIK_MSMILES)
        self.assertIsNotNone(metal_complex)
        # Ti + 3 Cl + bidentate iminopyridine + dimethylamine -> several ligands.
        self.assertGreaterEqual(len(metal_complex.ligands), 4)


if __name__ == "__main__":
    unittest.main()
