import unittest
import numpy as np
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles.utils.oin_aligner import OINDiscreteAligner

class TestOINAligner(unittest.TestCase):
    def test_linear_alignment(self):
        # Metal at origin
        metal_coords = np.array([0.0, 0.0, 0.0])
        
        # Ligand 1: Positive Z axis (LIN Slot 0: 0,0,1)
        # Ligand 2: Negative Z axis (LIN Slot 1: 0,0,-1)
        
        ligands = [
            {'smiles': '[Pt]', 'metal_coords': metal_coords}, # Metal is index 0
            {
                'smiles': '[Cl]', 
                'binding_atoms': [(1, 35.45, np.array([0.0, 0.0, 2.0]), 0)] 
            },
            {
                'smiles': '[Cl]', 
                'binding_atoms': [(2, 35.45, np.array([0.0, 0.0, -2.0]), 0)]
            }
        ]
        
        aligner = OINDiscreteAligner(metal_idx=0, ligands=ligands)
        oin_block = aligner.generate_canonical_vectors()
        
        # Expected: g:LIN|w:...
        self.assertIn("g:LIN", oin_block)
        # Check w tag. Should assign Slot 0 to Lig 1 and Slot 1 to Lig 2 (or vice versa aligned)
        self.assertIn("w:", oin_block)

    def test_square_planar_alignment(self):
        # Metal at origin
        metal_coords = np.array([0.0, 0.0, 0.0])
        
        # 4 Ligands in XY plane (SPL)
        # Lig 1: +X (Slot 0)
        # Lig 2: +Y (Slot 1)
        # Lig 3: -X (Slot 2)
        # Lig 4: -Y (Slot 3)
        
        ligands = [
            {'smiles': '[Pt]', 'metal_coords': metal_coords},
            {'smiles': 'L1', 'binding_atoms': [(1, 12.0, np.array([2.0, 0.0, 0.0]), 0)]},
            {'smiles': 'L2', 'binding_atoms': [(2, 12.0, np.array([0.0, 2.0, 0.0]), 0)]},
            {'smiles': 'L3', 'binding_atoms': [(3, 12.0, np.array([-2.0, 0.0, 0.0]), 0)]},
            {'smiles': 'L4', 'binding_atoms': [(4, 12.0, np.array([0.0, -2.0, 0.0]), 0)]},
        ]
        
        aligner = OINDiscreteAligner(metal_idx=0, ligands=ligands)
        oin_block = aligner.generate_canonical_vectors()
        
        self.assertIn("g:SPL", oin_block)

    def test_haptic_group(self):
        # 2 Ligands: 1 Monodentate, 1 Bidentate (treated as haptic/chelate?)
        # Aligner logic: "reduce hapticity" -> Groups binding atoms by distance.
        # Let's create a "bidentate" ligand with 2 atoms close to each other binding to metal?
        # Or 2 separate binding atoms in one ligand frag.
        
        metal_coords = np.array([0.0, 0.0, 0.0])
        
        # Lig A: +Z
        # Lig B: -Z, but has 2 atoms binding close to each other (Haptic?)
        
        ligands = [
            {'smiles': 'M', 'metal_coords': metal_coords},
            # Lig 1
            {'smiles': 'L1', 'binding_atoms': [(1, 14.0, np.array([0.0, 0.0, 2.0]), 0)]},
            # Lig 2 (Bidentate/Haptic)
            {
                'smiles': 'L2', 
                'binding_atoms': [
                    (2, 12.0, np.array([0.1, 0.0, -2.0]), 0),
                    (3, 12.0, np.array([-0.1, 0.0, -2.0]), 1)
                ]
            }
        ]
        
        # Distance between Lig2 atoms: 0.2 < 1.6 threshold. Should be grouped.
        # N_eff = 2. Should be LIN.
        
        aligner = OINDiscreteAligner(metal_idx=0, ligands=ligands)
        oin_block = aligner.generate_canonical_vectors()
        
        self.assertIn("g:LIN", oin_block)
        # Should contain references to atoms of Lig 2
        # w tag should have Rank 2, indices 0 and 1 assigned to same Slot.
        
        # Just check it runs and produces output
        self.assertIn("w:", oin_block)

if __name__ == '__main__':
    unittest.main()
