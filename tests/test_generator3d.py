import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from oinsmiles.generator3d import generate_3d_structures, get_xyz_string


def test_generation():
    msmiles = "[Zr+4]|[Cl-:2]|[Cl-:3]|[N:1]1=C(C[C-:4]2[CH:4]=[CH:4][CH:4]=[CH:4]2)C=CC=C1(C[C-:5]3[CH:5]=[CH:5][CH:5]=[CH:5]3)|5_trigonal_bipyramidal"
    print(f"Testing generation for:\n{msmiles}")

    mols = generate_3d_structures(msmiles, num_conformers=1)

    assert len(mols) > 0, "Failed to generate any conformers"

    print("Generation successful!")
    print("XYZ Format:")
    print(get_xyz_string(mols[0]))


if __name__ == "__main__":
    test_generation()
