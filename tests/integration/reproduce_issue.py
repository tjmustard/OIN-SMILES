import logging
import os
import sys
import tempfile

# Configure Logging
logging.basicConfig(level=logging.INFO)

# Setup path to include src
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../"))
src_path = os.path.join(project_root, "src")
sys.path.append(src_path)

from oinsmiles.utils.xyz2mol import get_oin_string, get_tmc_mol

trans_platin_xyz = """11
Transplatin
Pt     3.4456590331    0.5282513145    1.5566718340
Cl     5.3291166525   -0.4152532436    0.6270972055
Cl     1.5655270402    1.4634462578    2.5017817097
N      4.5660353118    2.0257535155    2.3862004479
H      3.9431541405    2.6312141638    2.9275292703
H      5.2834367779    1.6272738219    2.9951065752
N      2.3209838335   -0.9683094799    0.7307241558
H      2.9300914940   -1.7630556606    0.5207702333
H      1.5778049821   -1.2454955046    1.3751547164
H      1.8962264678   -0.6338889611   -0.1367039989
H      5.0271597487    2.5669933559    1.6522240323
"""

with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as tmp:
    tmp.write(trans_platin_xyz)
    tmp_path = tmp.name

try:
    print("Running conversion for TransPlatin...")
    mol, coords = get_tmc_mol(tmp_path, 0)
    oin = get_oin_string(mol, coords)
    print(f"OIN: {oin}")

    # Expected: Trans: [Pt_SPL].[Cl][0].N[1].[Cl][2].N[3]
    # Current behavior might be mass-sorted.
    expected = "[Pt_SPL].[Cl][0].N[1].[Cl][2].N[3]"

    if oin == expected:
        print("SUCCESS: Matches expected")
    else:
        print(f"FAILURE: Expected {expected}")
        print(f"         Got      {oin}")

finally:
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
