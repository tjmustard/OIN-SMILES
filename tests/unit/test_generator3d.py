import unittest

from oinsmiles.generator3d import generate_3d_structures, get_xyz_string


class TestGenerator3D(unittest.TestCase):
    def test_generation(self):
        msmiles = "[Zr+4]|[Cl-:2]|[Cl-:3]|[N:1]1=C(C[C-:4]2[CH:4]=[CH:4][CH:4]=[CH:4]2)C=CC=C1(C[C-:5]3[CH:5]=[CH:5][CH:5]=[CH:5]3)|5_trigonal_bipyramidal"

        mols = generate_3d_structures(msmiles, num_conformers=1)

        self.assertGreater(len(mols), 0, "Failed to generate any conformers")
        xyz = get_xyz_string(mols[0])
        self.assertIn("Zr", xyz)


if __name__ == "__main__":
    unittest.main()
