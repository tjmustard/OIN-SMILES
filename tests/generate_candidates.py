import os
import sys
from rdkit import Chem
from oinsmiles.generation.engine import OIN3DGenerator
from oinsmiles.core.translator import XYZToSMILES

def generate_cisplatin_candidate():
    # US-001/Test 6 of Molassembler MiniPRD
    oin = "[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}"
    gen = OIN3DGenerator()
    xyz = gen.generate(oin)
    
    os.makedirs("tests/candidate_outputs", exist_ok=True)
    with open("tests/candidate_outputs/molassembler_cisplatin.xyz", "w") as f:
        f.write(xyz)
    print("Generated tests/candidate_outputs/molassembler_cisplatin.xyz")

def generate_chiral_candidates():
    # For MiniPRD_ChiralTests Novel Test 5
    translator = XYZToSMILES()
    
    chiral_files = [
        "tests/integration/PdCl2-R-BINAP.xyz",
        "tests/integration/PdCl2-RR-BDNN.xyz",
        "tests/integration/PdCl2-RR-BDPP.xyz"
    ]
    
    for fpath in chiral_files:
        if os.path.exists(fpath):
            oin = translator.convert(fpath)
            # Save to candidate_outputs
            out_name = os.path.basename(fpath).replace(".xyz", "_encoded.smi")
            with open(f"tests/candidate_outputs/{out_name}", "w") as f:
                f.write(oin)
            print(f"Generated tests/candidate_outputs/{out_name}")
        else:
            print(f"Skipping {fpath} (not found)")

if __name__ == "__main__":
    generate_cisplatin_candidate()
    generate_chiral_candidates()
