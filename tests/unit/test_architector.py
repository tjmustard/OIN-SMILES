import unittest
from unittest.mock import MagicMock, patch
import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.generation.architector_adapter import ArchitectorAdapter
from oinsmiles.generation.oin_parser import ParsedOIN, OINVector

class TestArchitectorAdapter(unittest.TestCase):
    def setUp(self):
        self.adapter = ArchitectorAdapter()

    @patch('oinsmiles.generation.architector_adapter.Chem')
    @patch('oinsmiles.generation.architector_adapter.element') # Mock mendeleev
    def test_convert(self, mock_element, mock_chem):
        # Setup inputs
        # Use valid slot_idx for SPL (0..3)
        vec = OINVector(atom_idx=-1, vector=(1.0, 0.0, 0.0), fragment_idx=1, atom_in_fragment_idx=0, slot_idx=0)
        parsed = ParsedOIN(
            smiles="[Pt].[Cl]",
            fragments=["[Pt]", "[Cl]"],
            metal_fragment_idx=0,
            vectors=[vec],
            original_oin="raw",
            geometry="SPL"
        )

        # Setup RDKit mocks
        mock_mol_metal = MagicMock()
        mock_mol_metal.GetNumAtoms.return_value = 1
        mock_mol_metal.GetAtomWithIdx.return_value.GetSymbol.return_value = "Pt"
        
        mock_mol_ligand = MagicMock()
        mock_mol_ligand.GetNumAtoms.return_value = 1
        mock_atom = mock_mol_ligand.GetAtomWithIdx.return_value
        mock_atom.GetSymbol.return_value = "Cl"
        mock_atom.GetNumExplicitHs.return_value = 0
        mock_atom.GetNumImplicitHs.return_value = 0
        mock_atom.GetNoImplicit.return_value = False
        mock_atom.GetDegree.return_value = 0

        # Side effect for MolFromSmiles
        def mol_from_smiles(s, **kwargs):
            if "Pt" in s: return mock_mol_metal
            if "Cl" in s: return mock_mol_ligand
            return None
        mock_chem.MolFromSmiles.side_effect = mol_from_smiles
        
        # Mock MolToSmiles to return expected string
        mock_chem.MolToSmiles.return_value = "[Cl]"

        # Setup Mendeleev mocks
        mock_elem_obj = MagicMock()
        mock_elem_obj.covalent_radius_pyykko = 100.0 # 1.0 Angstrom
        mock_element.return_value = mock_elem_obj

        # Execute
        args = self.adapter.convert(parsed)

        # Verify
        # Verify
        self.assertEqual(args['core']['metal'], "Pt")
        self.assertEqual(len(args['ligands']), 1)
        
        lig = args['ligands'][0]
        self.assertEqual(lig['smiles'], "[Cl]")
        self.assertEqual(lig['coordinating_atoms'], [0])
        
        # Check site coords
        site_coords = args['core']['coordList']
        self.assertEqual(len(site_coords), 1)
        # Distance = (100+100)/100 = 2.0 Angstrom
        # Vector (1,0,0) -> Pos (2.0, 0.0, 0.0)
        self.assertAlmostEqual(site_coords[0][0], 2.0)

if __name__ == '__main__':
    unittest.main()
