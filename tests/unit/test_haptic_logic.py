import unittest
import numpy as np
from oinsmiles.generation.haptic import HapticTransformer, HapticResolver, TEMPLATES
from oinsmiles.generation.architector_adapter import ArchitectorAdapter
from oinsmiles.generation.oin_parser import ParsedOIN, OINVector

class TestHapticLogic(unittest.TestCase):
    
    def test_transformer_math_n5(self):
        """Test N=5 logic (Pentagon)"""
        z = np.array([0,0,1]) # Up
        x = np.array([1,0,0]) # Ref
        
        vectors = HapticTransformer.transform_vectors(5, z, x, cone_spread=0.2)
        self.assertEqual(len(vectors), 5)
        
        # V0 should be roughly Z + spread*X
        v0 = vectors[0]
        # Check angle with Z is roughly consistent
        dot_z = np.dot(v0, z)
        # expected angle? tan(alpha) = 0.2? approx.
        # normalize(Z + 0.2*X) -> (0,0,1) + (0.2,0,0) = (0.2, 0, 1). Norm = sqrt(1.04).
        # v0 = (0.2, 0, 1) / 1.02.
        self.assertTrue(v0[2] > 0.9)
        self.assertTrue(v0[0] > 0.1) # Positive X
        
        # V1 should be rotated 72 deg
        v1 = vectors[1]
        # In projection, angle between v0_xy and v1_xy should be 72 deg
        # v0_xy = (0.2, 0)
        # v1_xy should be (0.2 * cos72, 0.2 * sin72)
        v1_xy = v1[:2]
        angle = np.degrees(np.arctan2(v1[1], v1[0]))
        self.assertAlmostEqual(angle, 72.0, delta=1.0)

    def test_adapter_haptic_integration(self):
        """Verify Adapter uses Resolver correctly"""
        adapter = ArchitectorAdapter(scaling_factor=1.0)
        
        # Mock Input: [Fe].c1ccccc1{0} (Cp ring)
        # 5 vectors, all slot 0.
        # Let's say heading on index 0.
        
        vectors = []
        for i in range(5):
            vectors.append(OINVector(
                atom_idx=-1,
                vector=(0.0, 0.0, 1.0), # Dummy centroid
                fragment_idx=1, # Ligand
                atom_in_fragment_idx=i,
                haptic_heading=(i==2), # Atom 2 is Heading
                haptic_direction=-1, # Reverse
                slot_idx=0
            ))
            
        parsed = ParsedOIN(
            smiles="[Fe].c1ccccc1",
            fragments=["[Fe]", "c1ccccc1"],
            metal_fragment_idx=0,
            vectors=vectors,
            original_oin="...",
            geometry="LIN" # Use LIN for simplicity (Slot 0 is Z)
        )
        
        res = adapter.convert(parsed)
        coords = res['parameters']['site_coords']
        
        # Should have 5 coords
        self.assertEqual(len(coords), 5)
        
        # Verify they are DISTINCT
        # If logic failed, they would all be same (centroid)
        c0 = np.array(coords[0])
        c1 = np.array(coords[1])
        c2 = np.array(coords[2])
        
        dist = np.linalg.norm(c0 - c1)
        self.assertTrue(dist > 0.1, "Coordinates should be distinct")
        
        # Verify Heading Logic
        # Atom 2 was Heading. So c2 should be V0 (0 deg).
        # Atom 2 corresponds to index 2 in coords list (since sorted).
        # V0 for LIN slot 0 (Z, Ref X) is in XZ plane (y=0).
        self.assertAlmostEqual(c2[1], 0.0, delta=0.01) # Y should be 0
        self.assertTrue(c2[0] > 0) # Positive X
        
        # Atom 3: Delta = (3-2)*-1 = -1 -> Index 4 (V4 -> 288 deg -> -72 deg)
        # sin(-72) is negative Y.
        c3 = np.array(coords[3])
        self.assertTrue(c3[1] < 0) # Negative Y
        
        # Atom 1: Delta = (1-2)*-1 = 1 -> Index 1 (V1 -> 72 deg)
        # sin(72) is positive Y.
        c1 = np.array(coords[1])
        self.assertTrue(c1[1] > 0) # Positive Y

if __name__ == '__main__':
    unittest.main()
