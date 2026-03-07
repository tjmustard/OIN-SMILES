"""
Sanitization logic for OIN.
Locks chemical identity of ligands (Zone A) to prevent SMILES drift during geometry handling.
"""
from rdkit import Chem
from typing import List, Tuple
import re

def sanitize_metal_complex(mol: Chem.Mol, metal_idx: int = -1) -> Chem.Mol:
    """
    locks the explicit H count and disables implicit Hs for all atoms 
    bonded to the metal (Zone A) and the metal itself.
    
    Args:
        mol: The RDKit molecule (editable or not, we will make a writable copy/edit).
        metal_idx: Index of the metal atom. If -1, tries to find a Transition Metal.
        
    Returns:
        The sanitized molecule (modified in place if possible, but RDKit logic usually suggests editing RWMol).
    """
    # Create an editable molecule or work on the existing one if allowed.
    # Mol is usually passed as object. We'll modify it.
    
    if metal_idx == -1:
        # Simple heuristic: Find first transition metal? 
        # Or atomic number > 20 not in main group?
        # For OIN, we usually know the metal.
        # Let's search for the first metal atom.
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() >= 21: # Very rough check, but serves for now.
                metal_idx = atom.GetIdx()
                break
    
    if metal_idx == -1:
        # No metal found, nothing to sanitize in Zone A context
        return mol

    metal_atom = mol.GetAtomWithIdx(metal_idx)
    
    # 1. Identify Zone A (Metal + Neighbors)
    zone_a_indices = [metal_idx]
    for neighbor in metal_atom.GetNeighbors():
        zone_a_indices.append(neighbor.GetIdx())
        
    # 2. Lock Attributes
    for idx in zone_a_indices:
        atom = mol.GetAtomWithIdx(idx)
        
        # Lock H Count
        h_count = atom.GetTotalNumHs()
        atom.SetNumExplicitHs(h_count)
        atom.SetNoImplicit(True)
        
    # Standard sanitization to refresh properties (but keep explicit Hs)
    # Chem.SanitizeMol(mol) might strive to re-aromatize.
    # We want to PREVENT changes.
    mol.UpdatePropertyCache(strict=False)
    
    return mol

class OINSanitizer:
    @staticmethod
    def generate_robust_smiles(ligand_mol: Chem.Mol, binding_indices_in_ligand: List[int]) -> Tuple[str, Chem.Mol]:
        """
        Generates a canonical SMILES string where ALL binding atoms
        are forced to have explicit brackets. This prevents 'Drift'
        where c1 becomes [cH]1 or vice versa between runs.
        
        ligand_mol: RDKit Molecule Object of the ligand fragment
        binding_indices_in_ligand: List of atom indices (int) in the ligand 
                                   that bond to the metal.
        """
        # Create a modifiable copy
        rw_mol = Chem.RWMol(ligand_mol)
        
        # 1. Force Explicit H attributes on Zone A atoms
        for idx in binding_indices_in_ligand:
            atom = rw_mol.GetAtomWithIdx(idx)
            
            # Get current H count robustly to avoid Pre-condition Violations
            # If NoImplicit is set, we rely on Explicit count.
            if atom.GetNoImplicit():
                total_h = atom.GetNumExplicitHs()
            else:
                try:
                    total_h = atom.GetTotalNumHs()
                except RuntimeError:
                    # Fallback if property calculation fails
                    total_h = atom.GetNumExplicitHs()
            
            # Force this to be Explicit. 
            # This forces RDKit to write brackets like [cH] or [CH3] 
            # instead of c or C.
            atom.SetNumExplicitHs(total_h)
            atom.SetNoImplicit(True)
            
            # FORCE BRACKETS: Set a unique dummy atom map number.
            # RDKit ALWAYS brackets mapped atoms (e.g. [C:1]).
            # We use a unique ID (idx+1) to target specific atoms for post-processing cleanup.
            atom.SetAtomMapNum(idx + 1)
            
            # Update property cache to check valence state
            try:
                atom.GetOwningMol().UpdatePropertyCache(strict=False)
            except:
                pass
                
            # Check for valence deficit (to handle radicals like [Cl] vs Cl->HCl)
            # Only do this if not already set (though usually 0)
            if atom.GetNumRadicalElectrons() == 0:
                pt = Chem.GetPeriodicTable()
                default_val = pt.GetValenceList(atom.GetAtomicNum())[0]
                
                # Force Manual Calculation to avoid RDKit Pre-condition Violations
                # when properties are invalid/uncalculated
                explicit_h = float(atom.GetNumExplicitHs())
                bond_sum = 0.0
                for bond in atom.GetBonds():
                    bond_sum += bond.GetBondTypeAsDouble()
                
                current_val = explicit_h + bond_sum
                
                deficit = default_val - current_val
                
                if current_val < default_val:
                    # Ensure deficit is valid (non-negative int)
                    if deficit > 0:
                        atom.SetNumRadicalElectrons(int(deficit))
                        try:
                            atom.GetOwningMol().UpdatePropertyCache(strict=False)
                        except:
                            pass
            
        # 2. Restore Bond Stereo for Axial Chirality
        # These was captured in CIPAssigner.assign_all and stored on the bond.
        for bond in rw_mol.GetBonds():
            if bond.HasProp("OIN_BondStereo"):
                stereo_val = bond.GetProp("OIN_BondStereo")
                # Map back to Enum
                for member in Chem.BondStereo.values.values():
                    if str(member) == stereo_val:
                        bond.SetStereo(member)
                        break

        # 3. Generate Canonical SMILES
        # isomericSmiles=True ensures we keep stereochem info if present
        kmol = rw_mol.GetMol()
        try:
            kmol.UpdatePropertyCache(strict=False)
        except:
            pass
            
        smiles = Chem.MolToSmiles(kmol, isomericSmiles=True, canonical=True)
        
        # 3. Post-Processing Cleanup & Enforcement
        # RDKit can sometimes add 'H' to brackets (e.g. [cH:1]) even if ExplicitHs=0 was set,
        # especially for radicals or aromatic atoms.
        # We enforce strict "No H" if ExplicitHs was 0.
        
        for idx in binding_indices_in_ligand:
             atom = rw_mol.GetAtomWithIdx(idx)
             if atom.GetNumExplicitHs() == 0:
                 map_num = idx + 1
                 # Regex: Match [Symbol H : MapNum] and replace with [Symbol : MapNum]
                 # We capture the symbol ([a-zA-Z]+)
                 pattern = r'(\[[a-zA-Z]+)H:' + str(map_num) + r'\]'
                 replacement = r'\1:' + str(map_num) + r']'
                 
                 if re.search(pattern, smiles):
                     smiles = re.sub(pattern, replacement, smiles)

        # DO NOT REMOVE MAP NUMBERS YET
        # We need them for mapping back indices in xyz2mol
        # smiles = re.sub(r':\d+\]', ']', smiles)
        
        return smiles, kmol
