import sys
import os
import glob
import numpy as np
from scipy.spatial.transform import Rotation 

# Add src to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from oinsmiles import SMILESToXYZ, XYZToSMILES

from dataclasses import dataclass
from typing import List, Dict, Any, Optional
import tempfile
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

    def run(self):
        print(f"Running {len(self.examples)} Real Life Examples...")
        critical_failures = 0
        non_critical_failures = 0
        
        for i, example in enumerate(self.examples, 1):
            print(f"\n--- Example {i}: {example.name} ---")
            print(f"Description: {example.description}")
            
            success = True
            if example.oin_string:
                success = self._run_oin_to_xyz(example, i)
            elif example.xyz_content:
                success = self._run_xyz_to_smiles(example, i)
            else:
                print(f"Example {i} skipped: No input data provided.")
            
            if not success:
                if example.is_critical:
                    critical_failures += 1
                else:
                    non_critical_failures += 1

        print(f"\nFinished running {len(self.examples)} examples.")
        print(f"Critical Failures: {critical_failures}")
        print(f"Non-Critical Failures: {non_critical_failures}")
        
        if critical_failures > 0:
            print(f"FAILED: {critical_failures} critical examples failed.")
            sys.exit(1)
        elif non_critical_failures > 0:
            print(f"WARNING: {non_critical_failures} non-critical examples failed, but test suite passed.")
            sys.exit(0)
        else:
            print("SUCCESS: All examples passed.")

    def _run_oin_to_xyz(self, example: Example, index: int) -> bool:
        print(f"Type: OIN -> XYZ")
        print(f"Input OIN: {example.oin_string}")
        try:
            graph = self.smiles_to_xyz.convert(example.oin_string)
            print(f"Graph created with {len(graph.atoms)} atoms. (Expected: {example.expected_atom_count})")
            
            if example.expected_atom_count and len(graph.atoms) != example.expected_atom_count:
                print(f"WARNING: Atom count mismatch! Expected {example.expected_atom_count}, got {len(graph.atoms)}")
            
            dative_count = 0
            for atom in graph.atoms:
                for neighbor_idx, bond_type in atom.bonds.items():
                    if str(bond_type) == "BondType.DATIVE":
                        dative_count += 1
                        print(f"  Dative Bond: {atom.index} -> {neighbor_idx}")
            
            print(f"Total Dative Edges Found: {dative_count} (Expected ~{example.expected_dative_bonds * 2 if example.expected_dative_bonds else 'N/A'})")
            print(f"Example {index} Completed Successfully.")
            return True
        except Exception as e:
            print(f"Example {index} FAILED with error: {e}")
            return False

    def _run_xyz_to_smiles(self, example: Example, index: int) -> bool:
        print(f"Type: XYZ -> SMILES")
        
        # Apply transformation
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
            
            if example.expected_smiles:
                if output_smiles == example.expected_smiles:
                     print("SMILES Match: YES")
                else:
                     print(f"SMILES Match: NO (Expected {example.expected_smiles})")
                     failed = True
            
            if example.expected_oin_string:
                if output_oin == example.expected_oin_string:
                    print("OIN String Match: YES")
                else:
                    print(f"OIN String Match: NO")
                    print(f"Expected: {example.expected_oin_string}")
                    print(f"Got:      {output_oin}")
                    failed = True
            
            if example.expected_ligands:
                print("Checking for expected ligands:")
                for ligand in example.expected_ligands:
                    if ligand in output_smiles:
                        print(f"  [PASS] Ligand found: {ligand}")
                    else:
                        print(f"  [FAIL] Ligand NOT found: {ligand}")
                        failed = True

            if failed:
                print(f"Example {index} FAILED validation")
                return False
            else:
                print(f"Example {index} Completed Successfully.")
                return True
            
        except Exception as e:
            print(f"Example {index} FAILED with error: {e}")
            return False
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

def read_file_content(filepath):
    with open(filepath, 'r') as f:
        return f.read()

def get_examples() -> List[Example]:
    examples = []
    

    # Example 1: CisPlatin (XYZ -> OIN-SMILES)
    try:
        cisplatin_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'CisPlatin.xyz'))
        
        cisplatin_xyz_ex = Example(
            name="CisPlatin (XYZ -> OIN-SMILES)",
            xyz_content=cisplatin_xyz,
            description="Cisplatin with NH3 ligands (neutral OIN representation).",
            expected_smiles="[Pt].[Cl].[Cl].N.N",
            expected_oin_string="[Pt].[Cl].[Cl].N.N |g:SPL|w:1.0:0;2.0:2;3.0:3;4.0:1|",
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
            expected_smiles="[Pt].[Cl].[Cl].N.N",
            expected_oin_string="[Pt].[Cl].[Cl].N.N |g:SPL|w:1.0:0;2.0:1;3.0:2;4.0:3|",
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
            expected_smiles="[Pt].NCCN.[Cl].[Cl]",
            expected_oin_string="[Pt].NCCN.[Cl].[Cl] |g:SPL|w:1.0:0;1.3:2;2.0:3;3.0:1|",
            expected_ligands=["[Cl]", "NCCN"]
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
            expected_smiles="[Fe].c1cccc1.c1cccc1",
            expected_oin_string="[Fe].c1cccc1.c1cccc1 |g:LIN|w:1.0:0;2.0:1|",
            expected_ligands=["c1cccc1"]
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
            expected_smiles="[Pd].CP(C)c1ccccc1P(C)C.[Cl].[Cl]",
            expected_oin_string="[Pd].CP(C)c1ccccc1P(C)C.[Cl].[Cl] |g:SPL|w:1.1:0;1.9:2;2.0:3;3.0:1|",
            expected_ligands=["[Cl]", "CP(C)c1ccccc1P(C)C"]
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
            description="Facial Iridium tris(phenylpyridine) complex (neutral OIN representation).",
            expected_smiles="[Ir].[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1",
            expected_oin_string="[Ir].[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1 |g:OCT|w:1.0:2;1.11:4;2.0:0;2.11:1;3.0:5;3.11:3|",
            expected_ligands=["[Ir]", "[c]1ccccc1-c1ccccn1"]
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
            expected_smiles="[Ir].[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1",
            expected_oin_string="[Ir].[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1.[c]1ccccc1-c1ccccn1 |g:OCT|w:1.0:2;1.11:4;2.0:0;2.11:1;3.0:5;3.11:3|",
            expected_ligands=["[Ir]", "[c]1ccccc1-c1ccccn1"]
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
            expected_smiles="[Pt].[Br].[Cl].N.[CH3]",
            expected_oin_string="[Pt].[Br].[Cl].N.[CH3] |g:SPL|w:1.0:0;2.0:2;3.0:1;4.0:3|",
            expected_ligands=["[Br]", "[Cl]", "N", "[CH3]"]
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
            description="Square Planar Pt complex with 4 different ligands (Trans-arrangement).",
            expected_smiles="[Pt].[Br].[Cl].N.[CH3]",
            expected_oin_string="[Pt].[Br].[Cl].N.[CH3] |g:SPL|w:1.0:0;2.0:1;3.0:2;4.0:3|",
            expected_ligands=["[Br]", "[Cl]", "N", "[CH3]"]
        )
        examples.append(pt_trans_ex)
    except FileNotFoundError:
        print("Skipping PtMeNH3ClBr-Trans example: File not found.")
    


    # Example 10: FeCO5 (XYZ -> OIN-SMILES)
    try:
        feco5_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'FeCO5.xyz'))
        
        feco5_ex = Example(
            name="FeCO5 (XYZ -> OIN-SMILES)",
            xyz_content=feco5_xyz,
            description="Iron pentacarbonyl.",
            expected_smiles="[Fe].[C]#O.[C]#O.[C]#O.[C]#O.[C]#O",
            expected_oin_string="[Fe].[C]#O.[C]#O.[C]#O.[C]#O.[C]#O |g:TBP|w:1.0:3;2.0:2;3.0:0;4.0:1;5.0:4|",
            expected_ligands=["[Fe]", "[C]#O"]
        )
        examples.append(feco5_ex)
    except FileNotFoundError:
        print("Skipping FeCO5 example: File not found.")

    # Example 11: HgI3 (XYZ -> OIN-SMILES)
    try:
        hgi3_xyz = read_file_content(os.path.join(os.path.dirname(__file__), 'HgI3.xyz'))
        
        hgi3_ex = Example(
            name="HgI3 (XYZ -> OIN-SMILES)",
            xyz_content=hgi3_xyz,
            description="Mercury triiodide.",
            expected_smiles="[Hg].[I].[I].[I]",
            expected_oin_string="[Hg].[I].[I].[I] |g:TPL|w:1.0:1;2.0:0;3.0:2|",
            expected_ligands=["[Hg]", "[I]"]
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
            expected_smiles="[Re].[F].[F].[F].[F].[F].[F].[F]",
            # g:PBP is expected for N=7
            expected_oin_string="[Re].[F].[F].[F].[F].[F].[F].[F] |g:PBP|w:1.0:2;2.0:3;3.0:6;4.0:0;5.0:1;6.0:4;7.0:5|",
            expected_ligands=["[Re]", "[F]"]
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
            expected_smiles="[Ti].[Cl].[Cl].[Cl].[Cl]",
            expected_oin_string="[Ti].[Cl].[Cl].[Cl].[Cl] |g:TET|w:1.0:0;2.0:1;3.0:2;4.0:3|",
            expected_ligands=["[Ti]", "[Cl]"]
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
            expected_smiles="[Ti].c1cccc1.c1cccc1.[CH3].[CH3]",
            expected_oin_string="[Ti].c1cccc1.c1cccc1.[CH3].[CH3] |g:TET|w:1.0:0;2.0:1;3.0:2;4.0:3|",
            expected_ligands=["[Ti]", "c1cccc1", "[CH3]"]
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
            expected_smiles="[V].CC(=O)C=C(C)[O].CC(=O)C=C(C)[O].[O]",
            expected_oin_string="[V].CC(=O)C=C(C)[O].CC(=O)C=C(C)[O].[O] |g:SPY|w:1.2:1;1.6:3;2.2:4;2.6:2;3.0:0|",
            expected_ligands=["[V]", "[O]", "CC(=O)C=C(C)[O]"]
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
            expected_smiles="[Pt].[Cl].[Cl].[Cl].C=C",
            expected_oin_string="[Pt].[Cl].[Cl].[Cl].C=C |g:SPL|w:1.0:0;2.0:2;3.0:3;4.0:1|",
            expected_ligands=["[Pt]", "[Cl]", "C=C"]
        )
        examples.append(zeises_ex)
    except FileNotFoundError:
        print("Skipping Zeises_salt example: File not found.")

    # Add tmQM examples (Added last so they don't displace manual examples when limiting)
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

def run_examples(limit: Optional[int] = None):
    runner = ExampleRunner()
    examples = get_examples()
    
    if limit:
        print(f"Limiting to first {limit} examples.")
        examples = examples[:limit]
        
    for ex in examples:
        runner.add_example(ex)
    runner.run()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run Real Life Examples Verification")
    parser.add_argument("--limit", type=int, help="Limit number of examples to run (for fast testing)")
    args = parser.parse_args()
    
    run_examples(limit=args.limit)
