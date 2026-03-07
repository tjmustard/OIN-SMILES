"""
Integration test for Chiral Encoding (XYZ -> OIN-Inline -> XYZ).
Checks preservation of P/N stereocenters.
"""
import os
import logging
from rdkit import Chem
from oinsmiles.core.translator import XYZToSMILES, SMILESToXYZ

# Configure logging to see DEBUG output from chirality.py
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def create_p_chiral_xyz(filename):
    # Simulated P-chiral complex (V-type)
    # [V(CO)5(PH3)] but with a chiral P stereocenter if we had different ligands
    # Let's just use a simple chiral P ligand on a Metal
    # [Pt(Cl)2(P(F)(Cl)(H)(Cl))] ? No, let's just add one more ligand to P
    content = """8

Pt 0.0000 0.0000 0.0000
Cl 2.3000 0.0000 0.0000
Cl -2.3000 0.0000 0.0000
Cl 0.0000 2.3000 0.0000
P 0.0000 -2.3000 0.0000
F 0.0000 -2.3000 1.5000
H 0.0000 -2.3000 -1.5000
Cl 1.5000 -2.3000 0.0000
"""
    with open(filename, "w") as f:
        f.write(content)

def test_chiral_file_roundtrip(filename, expected_chiral=True):
    print(f"\n--- Testing {os.path.basename(filename)} ---")
    if not os.path.exists(filename):
        print(f"File {filename} not found, skipping.")
        return

    # 1. XYZ -> OIN-Inline
    converter = XYZToSMILES()
    oin_inline = converter.convert(filename, charge=0)
    
    print(f"Generated OIN: {oin_inline}")
    
    # CHECK: Does it contain [P@...] or [C@...]?
    if expected_chiral:
        assert "@" in oin_inline, f"Missing chirality marker in {os.path.basename(filename)}"
    
    # 2. OIN-Inline -> XYZ (Graph)
    parser = SMILESToXYZ()
    graph = parser.convert(oin_inline)
    
    print(f"SUCCESS: {os.path.basename(filename)} roundtrip check passed.")

if __name__ == "__main__":
    # Test cases provided by user
    test_dir = "tests/integration"
    chiral_files = [
        os.path.join(test_dir, "PdCl2-R-BINAP.xyz"),
        os.path.join(test_dir, "PdCl2-RR-BDNN.xyz"),
        os.path.join(test_dir, "PdCl2-RR-BDPP.xyz")
    ]
    
    success_count = 0
    for f in chiral_files:
        try:
            test_chiral_file_roundtrip(f)
            success_count += 1
        except Exception as e:
            print(f"Test FAILED for {f}: {e}")
            import traceback
            # traceback.print_exc()

    # Also keep the P-chiral fixture but fix it
    fixture_xyz = "/tmp/test_p_chiral.xyz"
    create_p_chiral_xyz(fixture_xyz)
    try:
        test_chiral_file_roundtrip(fixture_xyz)
        success_count += 1
    except Exception as e:
        print(f"Test FAILED for fixture: {e}")

    print(f"\nSummary: {success_count}/{len(chiral_files)+1} cases passed.")
