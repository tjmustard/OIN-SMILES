import sys
import os
import glob
import numpy as np
import traceback
from scipy.spatial.transform import Rotation 

# Add src to path
# INSERT BEFORE SYS.PATH APPEND TO OVERRIDE VENV? 
# No, usually we want to insert at 0.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

import oinsmiles
print(f"DEBUG: oinsmiles loaded from: {oinsmiles.__file__}")

from oinsmiles import SMILESToXYZ, XYZToSMILES

from dataclasses import dataclass
from typing import List, Dict, Any, Optional
from typing import List, Dict, Any, Optional
import tempfile
from reporting import VerificationReporter
try:
    from tmqm_expected import EXPECTED_TMQM_RESULTS
except ImportError:
    EXPECTED_TMQM_RESULTS = {}


@dataclass
class Example:
    name: str
    description: str
    # For OIN -> XYZ
    oin_string: Optional[str] = None
    expected_atom_count: Optional[int] = None
    expected_dative_bonds: Optional[int] = None
    # For XYZ -> SMILES
    xyz_content: Optional[str] = None
    expected_oin_string: Optional[str] = None
    expected_smiles: Optional[str] = None
    expected_ligands: Optional[List[str]] = None
    # Test control
    is_critical: bool = True
    fixed_orientation: bool = False

def transform_xyz_content(xyz_content: str) -> str:
    """
    Applies a random rotation and random translation (max +/- 10) to the XYZ content.
    """
    lines = xyz_content.strip().splitlines()
    if not lines:
        return xyz_content
        
    try:
        atom_count = int(lines[0].strip())
        comment = lines[1]
        atom_lines = lines[2:]
        
        # Parse atoms
        elements = []
        coords = []
        for line in atom_lines:
            parts = line.split()
            if not parts: continue
            elements.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
            
        coords_arr = np.array(coords)
        
        # 1. Random Rotation
        rot = Rotation.random()
        rotated_coords = rot.apply(coords_arr)
        
        # 2. Random Translation (+/- 10)
        translation = np.random.uniform(-10, 10, size=3)
        translated_coords = rotated_coords + translation
        
        # Reconstruct XYZ
        new_lines = [f"{atom_count}", f"{comment} [Transformed: Rotated + Translated {translation}]"]
        for elem, coord in zip(elements, translated_coords):
            new_lines.append(f"{elem} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
            
        return "\n".join(new_lines)
        
    except Exception as e:
        print(f"WARNING: Failed to transform XYZ content: {e}")
        return xyz_content


class ExampleRunner:
    def __init__(self):
        self.smiles_to_xyz = SMILESToXYZ()
        self.xyz_to_smiles = XYZToSMILES()
        self.examples: List[Example] = []

    def add_example(self, example: Example):
        self.examples.append(example)

    def run(self, output_dir: Optional[str] = None) -> None:
        print(f"Running {len(self.examples)} Real Life Examples...")
        reporter = VerificationReporter("XYZ <-> OIN Verification Report")

        for i, example in enumerate(self.examples, 1):
            print(f"\n--- Example {i}: {example.name} ---")
            print(f"Description: {example.description}")

            if example.oin_string:
                self._run_oin_to_xyz(example, i, reporter)
            elif example.xyz_content:
                self._run_xyz_to_smiles(example, i, reporter)
            else:
                print(f"Example {i} skipped: No input data provided.")
                reporter.log_failure(example.name, "No input data provided")

        reporter.print_summary()
        if output_dir:
            json_path = os.path.join(output_dir, "summary_integration.json")
            reporter.write_summary_json(json_path)
            print(f"JSON summary written to: {json_path}")

    def _run_oin_to_xyz(self, example: Example, index: int, reporter: VerificationReporter) -> bool:
        print(f"Type: OIN -> XYZ")
        print(f"Input OIN: {example.oin_string}")
        test_name = f"{example.name} (OIN->XYZ)"
        try:
            graph = self.smiles_to_xyz.convert(example.oin_string)
            print(f"Graph created with {len(graph.atoms)} atoms. (Expected: {example.expected_atom_count})")
            
            details = []
            failed = False
            
            if example.expected_atom_count and len(graph.atoms) != example.expected_atom_count:
                msg = f"Atom count mismatch! Expected {example.expected_atom_count}, got {len(graph.atoms)}"
                print(f"WARNING: {msg}")
                details.append(msg)
                failed = True
            
            dative_count = 0
            for atom in graph.atoms:
                for neighbor_idx, bond_type in atom.bonds.items():
                    if str(bond_type) == "BondType.DATIVE":
                        dative_count += 1
                        print(f"  Dative Bond: {atom.index} -> {neighbor_idx}")
            
            # Check dative bonds if specified? Logic missing in original but we can log success
            print(f"Total Dative Edges Found: {dative_count} (Expected ~{example.expected_dative_bonds * 2 if example.expected_dative_bonds else 'N/A'})")
            print(f"Example {index} Completed Successfully.")
            
            if failed:
                reporter.log_failure(test_name, "Validation Failed", got="\n".join(details))
                return False
            else:
                reporter.log_success(test_name, f"Atoms: {len(graph.atoms)}")
                return True
                
        except Exception as e:
            print(f"Example {index} FAILED with error: {e}")
            traceback.print_exc()
            reporter.log_failure(test_name, f"Exception: {str(e)}")
            return False

    def _run_xyz_to_smiles(self, example: Example, index: int, reporter: VerificationReporter) -> bool:
        print(f"Type: XYZ -> SMILES")
        test_name = f"{example.name} (XYZ->SMILES)"
        
        # Apply transformation
        if example.fixed_orientation:
            print("Skipping random rotation (Fixed Orientation requested)")
            transformed_content = example.xyz_content
        else:
            transformed_content = transform_xyz_content(example.xyz_content)
        
        # Create temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as tmp:
            tmp.write(transformed_content)
            tmp_path = tmp.name
        
        print(f"Created temp XYZ file: {tmp_path}")
        try:
            # Print first few lines of XYZ for verification
            print("XYZ Content Preview:")
            print("\n".join(example.xyz_content.splitlines()[:5]))
            
            output_oin = self.xyz_to_smiles.convert(tmp_path)
            print(f"Output OIN: {output_oin}")
            
            output_smiles = output_oin.split("|")[0].strip()
            print(f"Extracted SMILES: {output_smiles}")

            failed = False
            fail_reason = ""
            expected_val = ""
            got_val = ""
            
            if example.expected_smiles:
                if output_smiles != example.expected_smiles.strip():
                     print(f"SMILES Match: NO (Expected {example.expected_smiles})")
                     failed = True
                     fail_reason = "SMILES Mismatch"
                     expected_val = example.expected_smiles
                     got_val = output_smiles
            
            if example.expected_oin_string:
                if output_oin.strip() != example.expected_oin_string.strip():
                    print(f"OIN String Match: NO")
                    print(f"Expected: {example.expected_oin_string}")
                    print(f"Got:      {output_oin}")
                    if not failed: # Capture first failure
                        failed = True
                        fail_reason = "OIN String Mismatch"
                        expected_val = example.expected_oin_string
                        got_val = output_oin
            
            if example.expected_ligands:
                print("Checking for expected ligands:")
                missing = []
                for ligand in example.expected_ligands:
                    if ligand not in output_smiles:
                        print(f"  [FAIL] Ligand NOT found: {ligand}")
                        missing.append(ligand)
                    else:
                        print(f"  [PASS] Ligand found: {ligand}")
                
                if missing:
                    failed = True
                    if not fail_reason:
                        fail_reason = f"Missing Ligands: {missing}"

            if failed:
                print(f"Example {index} FAILED validation")
                reporter.log_failure(test_name, fail_reason, expected=expected_val, got=got_val)
                return False
            else:
                print(f"Example {index} Completed Successfully.")
                reporter.log_success(test_name, f"OIN: {output_oin}")
                return True
            
        except Exception as e:
            print(f"Example {index} FAILED with error: {e}")
            traceback.print_exc()
            reporter.log_failure(test_name, f"Exception: {str(e)}")
            return False
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

def read_file_content(filepath):
    with open(filepath, 'r') as f:
        return f.read()

def get_examples(include_tmqm: bool = False) -> List[Example]:
    examples = []


    # Example 1: CisPlatin (XYZ -> OIN-SMILES)
    try:
        cisplatin_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'CisPlatin.xyz'))
        
        cisplatin_xyz_ex = Example(
            name="CisPlatin (XYZ -> OIN-SMILES)",
            xyz_content=cisplatin_xyz,
            description="Cisplatin with NH3 ligands (neutral OIN representation).",
            expected_smiles="[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
            expected_oin_string="[Pt_SPL].[Cl]{0}.[Cl]{1}.N{2}.N{3}",
            expected_ligands=["[Cl]", "N"]
        )
        examples.append(cisplatin_xyz_ex)
    except FileNotFoundError:
        print("Skipping CisPlatin example: File not found.")

    # Example 2: TransPlatin (XYZ -> OIN-SMILES)
    try:
        transplatin_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'TransPlatin.xyz'))
        
        transplatin_xyz_ex = Example(
            name="TransPlatin (XYZ -> OIN-SMILES)",
            xyz_content=transplatin_xyz,
            description="Transplatin with NH3 ligands (neutral OIN representation).",
            expected_smiles="[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}",
            expected_oin_string="[Pt_SPL].[Cl]{0}.N{1}.[Cl]{2}.N{3}",
            expected_ligands=["[Cl]", "N"]
        )
        examples.append(transplatin_xyz_ex)
    except FileNotFoundError:
        print("Skipping TransPlatin example: File not found.")

    # Example 3: Cis-PtCl2(en) (XYZ -> OIN-SMILES)
    try:
        cisptcl2en_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'Cis-PtCl2(en).xyz'))
        
        cisptcl2en_xyz_ex = Example(
            name="Cis-PtCl2(en) (XYZ -> OIN-SMILES)",
            xyz_content=cisptcl2en_xyz,
            description="Converting XYZ coordinates of Cisplatin to OIN-SMILES (neutral components).",
            expected_smiles="[Pt_SPL].[NH2]{0}CC[NH2]{1}.[Cl]{2}.[Cl]{3}",
            expected_oin_string="[Pt_SPL].[NH2]{0}CC[NH2]{1}.[Cl]{2}.[Cl]{3}",
            expected_ligands=["Cl", "NH2"]
        )
        examples.append(cisptcl2en_xyz_ex)
    except FileNotFoundError:
        print("Skipping Cis-PtCl2(en) example: File not found.")
    
    # Example 4: Ferrocene (XYZ -> OIN-SMILES)
    try:
        ferrocene_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'Ferrocene.xyz'))
        
        ferrocene_ex = Example(
            name="Ferrocene (XYZ -> OIN-SMILES)",
            xyz_content=ferrocene_xyz,
            description="Ferrocene with eclipsed Cp rings (neutral OIN representation).",
            expected_smiles="[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1",
            expected_oin_string="[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1",
            fixed_orientation=True
        )
        examples.append(ferrocene_ex)
    except FileNotFoundError:
        print("Skipping Ferrocene example: File not found.")
    
    # Example 5: PdCl2PhenPhosMe (XYZ -> OIN-SMILES)
    try:
        pd_phenphos_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'PdCl2PhenPhosMe.xyz'))
        
        pd_phenphos_ex = Example(
            name="PdCl2PhenPhosMe (XYZ -> OIN-SMILES)",
            xyz_content=pd_phenphos_xyz,
            description="Palladium complex with Cl and PhenPhosMe ligands.",
            expected_smiles="[Pd_SPL].CP{0}(C)c1ccccc1P{1}(C)C.[Cl]{2}.[Cl]{3}",
            expected_oin_string="[Pd_SPL].CP{0}(C)c1ccccc1P{1}(C)C.[Cl]{2}.[Cl]{3}"
        )
        examples.append(pd_phenphos_ex)
    except FileNotFoundError:
        print("Skipping PdCl2PhenPhosMe example: File not found.")

    # Example 6: fac-Ir(ppy)3 (XYZ -> OIN-SMILES)
    try:
        fac_ir_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'fac-Ir(ppy)3.xyz'))
        
        fac_ir_ex = Example(
            name="fac-Ir(ppy)3 (XYZ -> OIN-SMILES)",
            xyz_content=fac_ir_xyz,
            description="fac-Ir(ppy)3 with corrected OIN.",
            expected_smiles="[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1",
            expected_oin_string="[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{5}1ccccc1-c1ccccn{1}1.c{2}1ccccc1-c1ccccn{4}1"
        )
        examples.append(fac_ir_ex)
    except FileNotFoundError:
        print("Skipping fac-Ir(ppy)3 example: File not found.")

    # Example 7: mer-Ir(ppy)3 (XYZ -> OIN-SMILES)
    try:
        mer_ir_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'mer-Ir(ppy)3.xyz'))
        
        mer_ir_ex = Example(
            name="mer-Ir(ppy)3 (XYZ -> OIN-SMILES)",
            xyz_content=mer_ir_xyz,
            description="Meridional Iridium tris(phenylpyridine) complex (neutral OIN representation).",
            expected_smiles="[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1",
            expected_oin_string="[Ir_OCT].c{0}1ccccc1-c1ccccn{3}1.c{1}1ccccc1-c1ccccn{5}1.c{2}1ccccc1-c1ccccn{4}1"
        )
        examples.append(mer_ir_ex)
    except FileNotFoundError:
        print("Skipping mer-Ir(ppy)3 example: File not found.")

    # Example 8: PtMeNH3ClBr-Cis (XYZ -> OIN-SMILES)
    try:
        pt_cis_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'PtMeNH3ClBr-Cis.xyz'))
        
        pt_cis_ex = Example(
            name="PtMeNH3ClBr-Cis (XYZ -> OIN-SMILES)",
            xyz_content=pt_cis_xyz,
            description="Square Planar Pt complex with 4 different ligands (Cis-arrangement).",
            expected_smiles="[Pt_SPL].[Br]{0}.[Cl]{1}.N{2}.[CH3]{3}",
            expected_oin_string="[Pt_SPL].[Br]{0}.[Cl]{1}.N{2}.[CH3]{3}"
        )
        examples.append(pt_cis_ex)
    except FileNotFoundError:
        print("Skipping PtMeNH3ClBr-Cis example: File not found.")

    # Example 9: PtMeNH3ClBr-Trans (XYZ -> OIN-SMILES)
    try:
        pt_trans_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'PtMeNH3ClBr-Trans.xyz'))
        
        pt_trans_ex = Example(
            name="PtMeNH3ClBr-Trans (XYZ -> OIN-SMILES)",
            xyz_content=pt_trans_xyz,
            description="Mixed ligand Platinum complex.",
            expected_smiles="[Pt_SPL].[Br]{0}.N{1}.[Cl]{2}.[CH3]{3}",
            expected_oin_string="[Pt_SPL].[Br]{0}.N{1}.[Cl]{2}.[CH3]{3}"
        )
        examples.append(pt_trans_ex)
    except FileNotFoundError:
        print("Skipping PtMeNH3ClBr-Trans example: File not found.")
    
    # Example 10: CuCN2 (XYZ -> OIN-SMILES)
    try:
        cucn2_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'CuCN2.xyz'))
        
        cucn2_ex = Example(
            name="CuCN2 (XYZ -> OIN-SMILES)",
            xyz_content=cucn2_xyz,
            description="Linear Cu complex with 2 cyanide ligands.",
            expected_smiles="[Cu_LIN].C{0}#N.C{1}#N",
            expected_oin_string="[Cu_LIN].C{0}#N.C{1}#N"
        )
        examples.append(cucn2_ex)
    except FileNotFoundError:
        print("Skipping CuCN2 example: File not found.")


    # Example 11: FeCO5 (XYZ -> OIN-SMILES)
    try:
        feco5_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'FeCO5.xyz'))
        
        feco5_ex = Example(
            name="FeCO5 (XYZ -> OIN-SMILES)",
            xyz_content=feco5_xyz,
            description="Iron pentacarbonyl.",
            expected_smiles="[Fe_TBP].C{0}#O.C{1}#O.C{2}#O.C{3}#O.C{4}#O",
            expected_oin_string="[Fe_TBP].C{0}#O.C{1}#O.C{2}#O.C{3}#O.C{4}#O"
        )
        examples.append(feco5_ex)
    except FileNotFoundError:
        print("Skipping FeCO5 example: File not found.")

    # Example 11: FeH2(CO)4 (XYZ -> OIN-SMILES)
    try:
        feh2co4_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'FeH2(CO)4.xyz'))

        feh2co4_ex = Example(
            name="FeH2(CO)4 (XYZ -> OIN-SMILES)",
            xyz_content=feh2co4_xyz,
            description="Iron dihydride tetracarbonyl.",
            expected_smiles="[Fe_OCT].[H]{0}.[H]{1}.C{2}#O.C{3}#O.C{4}#O.C{5}#O",
            expected_oin_string="[Fe_OCT].[H]{0}.[H]{1}.C{2}#O.C{3}#O.C{4}#O.C{5}#O"
        )
        examples.append(feh2co4_ex)
    except FileNotFoundError:
        print("Skipping FeH2(CO)4 example: File not found.")

    # Example 11: HgI3 (XYZ -> OIN-SMILES)
    try:
        hgi3_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'HgI3.xyz'))
        
        hgi3_ex = Example(
            name="HgI3 (XYZ -> OIN-SMILES)",
            xyz_content=hgi3_xyz,
            description="Mercury triiodide.",
            expected_smiles="[Hg_TPL].[I]{0}.[I]{1}.[I]{2}",
            expected_oin_string="[Hg_TPL].[I]{0}.[I]{1}.[I]{2}"
        )
        examples.append(hgi3_ex)
    except FileNotFoundError:
        print("Skipping HgI3 example: File not found.")

    # Example 12: ReF7 (XYZ -> OIN-SMILES)
    try:
        ref7_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'ReF7.xyz'))
        
        ref7_ex = Example(
            name="ReF7 (XYZ -> OIN-SMILES)",
            xyz_content=ref7_xyz,
            description="Rhenium heptafluoride.",
            expected_smiles="[Re_PBP].[F]{0}.[F]{1}.[F]{2}.[F]{3}.[F]{4}.[F]{5}.[F]{6}",
            expected_oin_string="[Re_PBP].[F]{0}.[F]{1}.[F]{2}.[F]{3}.[F]{4}.[F]{5}.[F]{6}"
        )
        examples.append(ref7_ex)
    except FileNotFoundError:
        print("Skipping ReF7 example: File not found.")

    # Example 13: TiCl4 (XYZ -> OIN-SMILES)
    try:
        ticl4_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'TiCl4.xyz'))
        
        ticl4_ex = Example(
            name="TiCl4 (XYZ -> OIN-SMILES)",
            xyz_content=ticl4_xyz,
            description="Titanium tetrachloride.",
            expected_smiles="[Ti_TET].[Cl]{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}",
            expected_oin_string="[Ti_TET].[Cl]{0}.[Cl]{1}.[Cl]{2}.[Cl]{3}"
        )
        examples.append(ticl4_ex)
    except FileNotFoundError:
        print("Skipping TiCl4 example: File not found.")

    # Example 14: TiCp2Me2 (XYZ -> OIN-SMILES)
    try:
        ticp2me2_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'TiCp2Me2.xyz'))
        
        ticp2me2_ex = Example(
            name="TiCp2Me2 (XYZ -> OIN-SMILES)",
            xyz_content=ticp2me2_xyz,
            description="Titanocene dimethyl.",
            expected_smiles="[Ti_TET].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}",
            expected_oin_string="[Ti_TET].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1<}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}",
            fixed_orientation=True
        )
        examples.append(ticp2me2_ex)
    except FileNotFoundError:
        print("Skipping TiCp2Me2 example: File not found.")

    # Example 15: VOacac2 (XYZ -> OIN-SMILES)
    try:
        voacac2_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'VOacac2.xyz'))
        
        # Note: acac appears as anionic radical form in OIN due to sanitization
        voacac2_ex = Example(
            name="VOacac2 (XYZ -> OIN-SMILES)",
            xyz_content=voacac2_xyz,
            description="Vanadyl acetylacetonate.",
            expected_smiles="[V_SPY].O{0}.CC(=O{1})C=C(C)O{4}.CC(=O{2})C=C(C)O{3}",
            expected_oin_string="[V_SPY].O{0}.CC(=O{1})C=C(C)O{4}.CC(=O{2})C=C(C)O{3}"
        )
        examples.append(voacac2_ex)
    except FileNotFoundError:
        print("Skipping VOacac2 example: File not found.")

    # Example 16: Zeises_salt (XYZ -> OIN-SMILES)
    try:
        zeises_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'Zeises_salt.xyz'))
        
        zeises_ex = Example(
            name="Zeises_salt (XYZ -> OIN-SMILES)",
            xyz_content=zeises_xyz,
            description="Zeise's Salt (PtCl3(C2H4)).",
            expected_smiles="[Pt_SPL].[Cl]{0}.[Cl]{1}.[CH2]{2>}=[CH2]{2}.[Cl]{3}",
            expected_oin_string="[Pt_SPL].[Cl]{0}.[Cl]{1}.[CH2]{2>}=[CH2]{2}.[Cl]{3}"
        )
        examples.append(zeises_ex)
    except FileNotFoundError:
        print("Skipping Zeises_salt example: File not found.")

    # TiCat1
    try:
        with open(os.path.join(os.path.dirname(__file__), 'TiCat1.xyz'), "r") as f:
            ticat1_xyz = f.read()
        ticat1_ex = Example( # Changed from XYZToOINExample to Example
            name="TiCat1 (XYZ -> OIN-SMILES)",
            xyz_content=ticat1_xyz,
            description="Titanium Catalyst 1",
            expected_smiles="[Ti_TET].C[Si](C)(c{0}1[cH]{0}[cH]{0}[cH]{0}[cH]{0<}1)c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}",
            expected_oin_string="[Ti_TET].C[Si](C)(c{0}1[cH]{0}[cH]{0}[cH]{0}[cH]{0<}1)c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}",
            fixed_orientation=True
        )
        examples.append(ticat1_ex)
    except FileNotFoundError:
        print("Skipping TiCat1 example: File not found.")

    # TiCat2
    try:
        with open(os.path.join(os.path.dirname(__file__), 'TiCat2.xyz'), "r") as f:
            ticat2_xyz = f.read()
        ticat2_ex = Example( # Changed from XYZToOINExample to Example
            name="TiCat2 (XYZ -> OIN-SMILES)",
            xyz_content=ticat2_xyz,
            description="Titanium Catalyst 2",
            expected_smiles="[Ti_TET].CC(C)(C)N{0}[Si](C)(C)c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}",
            expected_oin_string="[Ti_TET].CC(C)(C)N{0}[Si](C)(C)c{1}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1.[CH3]{2}.[CH3]{3}"
        )
        examples.append(ticat2_ex)
    except FileNotFoundError:
        print("Skipping TiCat2 example: File not found.")

    # TiCat3
    try:
        with open(os.path.join(os.path.dirname(__file__), 'TiCat3.xyz'), "r") as f:
            ticat3_xyz = f.read()
        ticat3_ex = Example(
            name="TiCat3 (XYZ -> OIN-SMILES)",
            xyz_content=ticat3_xyz,
            description="Titanium Catalyst 3",
            expected_smiles="[Ti_TPY].[CH3]{0}.[CH3]{1}.C[Si](C)(c{2}1[cH]{2}[cH]{2}c{2}2ccccc{2<}12)c{3}1[cH]{3}[cH]{3}c{3}2ccccc{3}12",
            expected_oin_string="[Ti_TPY].[CH3]{0}.[CH3]{1}.C[Si](C)(c{2}1[cH]{2}[cH]{2}c{2}2ccccc{2<}12)c{3}1[cH]{3}[cH]{3}c{3}2ccccc{3}12"
        )
        examples.append(ticat3_ex)
    except FileNotFoundError:
        print("Skipping TiCat3 example: File not found.")

    # TiCat4
    try:
        with open(os.path.join(os.path.dirname(__file__), 'TiCat4.xyz'), "r") as f:
            ticat4_xyz = f.read()
        ticat4_ex = Example(
            name="TiCat4 (XYZ -> OIN-SMILES)",
            xyz_content=ticat4_xyz,
            description="Titanium Catalyst 4",
            expected_smiles="[Ti_TPY].[CH3]{0}.[CH3]{1}.C[Si](C)(c{2}1[cH]{2>}[cH]{2}c{2}2ccccc{2}12)c{3}1[cH]{3}[cH]{3}c{3}2ccccc{3}12",
            expected_oin_string="[Ti_TPY].[CH3]{0}.[CH3]{1}.C[Si](C)(c{2}1[cH]{2>}[cH]{2}c{2}2ccccc{2}12)c{3}1[cH]{3}[cH]{3}c{3}2ccccc{3}12"
        )
        examples.append(ticat4_ex)
    except FileNotFoundError:
        print("Skipping TiCat4 example: File not found.")

    # Example: PdCl2-R-BINAP (axial-chiral BINAP ligand)
    try:
        pdcl2_binap_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'PdCl2-R-BINAP.xyz'))

        pdcl2_binap_ex = Example(
            name="PdCl2-R-BINAP (XYZ -> OIN-SMILES)",
            xyz_content=pdcl2_binap_xyz,
            description="Palladium complex with axial-chiral R-BINAP diphosphine ligand.",
            expected_smiles="[Pd_SPL].c1ccc(P{0}(c2ccccc2)c2ccc3ccccc3c2-c2c(P{1}(c3ccccc3)c3ccccc3)ccc3ccccc23)cc1.[Cl]{2}.[Cl]{3}",
            expected_oin_string="[Pd_SPL].c1ccc(P{0}(c2ccccc2)c2ccc3ccccc3c2-c2c(P{1}(c3ccccc3)c3ccccc3)ccc3ccccc23)cc1.[Cl]{2}.[Cl]{3}",
            expected_ligands=["P", "Cl"]
        )
        examples.append(pdcl2_binap_ex)
    except FileNotFoundError:
        print("Skipping PdCl2-R-BINAP example: File not found.")

    # Example: PdCl2-RR-BDNN (N-chiral diphosphine ligand)
    try:
        pdcl2_bdnn_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'PdCl2-RR-BDNN.xyz'))

        pdcl2_bdnn_ex = Example(
            name="PdCl2-RR-BDNN (XYZ -> OIN-SMILES)",
            xyz_content=pdcl2_bdnn_xyz,
            description="Palladium complex with RR-BDNN N-chiral diphosphine ligand (R,R-bis(diethylamino)naphthalene).",
            expected_smiles="[Pd_SPL].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}",
            expected_oin_string="[Pd_SPL].C[C@@H](C[C@H](C)N{0}(c1ccccc1)c1ccccc1)N{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}",
            expected_ligands=["N", "Cl"]
        )
        examples.append(pdcl2_bdnn_ex)
    except FileNotFoundError:
        print("Skipping PdCl2-RR-BDNN example: File not found.")

    # Example: PdCl2-RR-BDPP (P-chiral diphosphine ligand)
    try:
        pdcl2_bdpp_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'PdCl2-RR-BDPP.xyz'))

        pdcl2_bdpp_ex = Example(
            name="PdCl2-RR-BDPP (XYZ -> OIN-SMILES)",
            xyz_content=pdcl2_bdpp_xyz,
            description="Palladium complex with RR-BDPP P-chiral diphosphine ligand (R,R-bis(diphenylphosphinoyl)benzene).",
            expected_smiles="[Pd_SPL].C[C@@H](C[C@H](C)P{0}(c1ccccc1)c1ccccc1)P{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}",
            expected_oin_string="[Pd_SPL].C[C@@H](C[C@H](C)P{0}(c1ccccc1)c1ccccc1)P{1}(c1ccccc1)c1ccccc1.[Cl]{2}.[Cl]{3}",
            expected_ligands=["P", "Cl"]
        )
        examples.append(pdcl2_bdpp_ex)
    except FileNotFoundError:
        print("Skipping PdCl2-RR-BDPP example: File not found.")

    # Add tmQM examples (Added last so they don't displace manual examples when limiting)
    if include_tmqm:
        tmqm_dir = os.path.join(os.path.dirname(__file__), 'tmQM')
        if os.path.exists(tmqm_dir):
            tmqm_files = glob.glob(os.path.join(tmqm_dir, "*.xyz"))
            print(f"Found {len(tmqm_files)} tmQM examples.")
            for xyz_file in tmqm_files:
                try:
                    content = read_file_content(xyz_file)
                    filename = os.path.basename(xyz_file)
                    expected_data = EXPECTED_TMQM_RESULTS.get(filename)

                    tmqm_ex = Example(
                        name=f"tmQM - {filename}",
                        xyz_content=content,
                        description=f"tmQM example from {filename}",
                        is_critical=True if expected_data else False,
                        expected_oin_string=expected_data['oin'] if expected_data else None,
                        expected_smiles=expected_data['smiles'] if expected_data else None
                    )
                    examples.append(tmqm_ex)
                except Exception as e:
                    print(f"Error reading {xyz_file}: {e}")
        else:
            print(f"tmQM directory not found at {tmqm_dir}")

    return examples

def run_examples(limit: Optional[int] = None, output_dir: Optional[str] = None, include_tmqm: bool = False) -> None:
    runner = ExampleRunner()
    examples = get_examples(include_tmqm=include_tmqm)

    if limit:
        print(f"Limiting to first {limit} examples.")
        examples = examples[:limit]

    for ex in examples:
        runner.add_example(ex)
    runner.run(output_dir=output_dir)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run Real Life Examples Verification")
    parser.add_argument("--limit", type=int, help="Limit number of examples to run (for fast testing)")
    parser.add_argument("--output-dir", type=str, help="Directory to save JSON summary artifact")
    parser.add_argument("--include-tmqm", action="store_true", help="Include tmQM examples (slow)")
    args = parser.parse_args()

    run_examples(limit=args.limit, output_dir=args.output_dir, include_tmqm=args.include_tmqm)
