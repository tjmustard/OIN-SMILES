"""
Core logic for preserving and recovering chirality (P/N stereocenters and axial chirality) 
through fragmentation and OIN serialization.
"""
import logging
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
from typing import List, Optional, Tuple, Any, Dict

logger = logging.getLogger(__name__)

class CIPAssigner:
    """
    Handles robust CIP assignment and property capture on the connected molecule.
    Used BEFORE fragmentation.
    """
    @staticmethod
    def assign_all(mol: Chem.Mol, coords: Optional[List[List[float]]] = None) -> None:
        """
        Calculates CIP codes and stores them as persistent properties on atoms.
        Uses 3D coordinates if provided for perception.
        
        Propagates Sanitization exceptions to caller per MiniPRD US-005.
        """
        if mol is None:
            raise ValueError("Molecule object cannot be None for CIP assignment.")

        # 1. Add Conformer if coords are provided
        if coords:
            try:
                conf = Chem.Conformer(mol.GetNumAtoms())
                for i, pos in enumerate(coords):
                    if i < mol.GetNumAtoms():
                        conf.SetAtomPosition(i, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))
                mol.RemoveAllConformers()
                mol.AddConformer(conf)
            except Exception as e:
                logger.warning(f"Failed to add 3D conformer for stereochem perception: {e}")

        # 2. Extract coords from conformer if not provided
        if not coords and mol.GetNumConformers() > 0:
            conf = mol.GetConformer()
            coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]

        # 3. Ensure properties are calculated
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        
        # 4. Handle 3-coordinate P/N (Lone Pair Pseudo-atom Strategy)
        # RDKit often ignores 3rd-coordinate inversion centers if they lack a 4th neighbor.
        added_pseudo_indices = []
        if isinstance(mol, Chem.RWMol) and mol.GetNumConformers() > 0:
            added_pseudo_indices = PseudoAtomStrategy.add_lone_pair_pseudos(mol)

        # 5. Perceive status from 3D
        if mol.GetNumConformers() > 0:
            Chem.AssignStereochemistryFrom3D(mol)
        
        # We use force=True to ensure we don't rely on stale tags
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

        # 6. Capture _CIPCode onto atoms as persistent string props
        for atom in mol.GetAtoms():
            if atom.HasProp("_CIPCode"):
                cip = atom.GetProp("_CIPCode")
                atom.SetProp("OIN_CIP", cip)
                logger.debug(f"Captured CIP {cip} for atom {atom.GetSymbol()} at idx {atom.GetIdx()}")

        # 7. Manual Axial Chirality Perception (Biaryls)
        # Standard FindPotentialStereoBonds often misses single-bond atropisomers
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing():
                u, v = bond.GetBeginAtom(), bond.GetEndAtom()
                if u.GetIsAromatic() and v.GetIsAromatic():
                    if coords is not None:
                         u_neighbors = [n.GetIdx() for n in u.GetNeighbors() if n.GetIdx() != v.GetIdx()]
                         v_neighbors = [n.GetIdx() for n in v.GetNeighbors() if n.GetIdx() != u.GetIdx()]
                         if u_neighbors and v_neighbors:
                             idx1, idx2, idx3, idx4 = u_neighbors[0], u.GetIdx(), v.GetIdx(), v_neighbors[0]
                             p1, p2, p3, p4 = [np.array(coords[i], dtype=float) for i in [idx1, idx2, idx3, idx4]]
                             
                             b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
                             n1 = np.cross(b1, b2).astype(float)
                             n2 = np.cross(b2, b3).astype(float)
                             norm1, norm2 = np.linalg.norm(n1), np.linalg.norm(n2)
                             if norm1 > 1e-6: n1 /= norm1
                             if norm2 > 1e-6: n2 /= norm2
                             
                             norm_b2 = np.linalg.norm(b2)
                             if norm_b2 > 1e-6:
                                 m1 = np.cross(n1, b2 / norm_b2)
                                 x, y = np.dot(n1, n2), np.dot(m1, n2)
                                 dihedral = np.arctan2(y, x)
                                 
                                 stereo_str = "STEREOATROP_CW" if dihedral > 0 else "STEREOATROP_CCW"
                                 bond.SetProp("OIN_BondStereo", stereo_str)
                                 logger.debug(f"Detected axial chirality at bond {u.GetIdx()}-{v.GetIdx()}: {stereo_str}")

        # Standard double bond perception
        Chem.FindPotentialStereoBonds(mol)
        for bond in mol.GetBonds():
            if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                 if not bond.HasProp("OIN_BondStereo"):
                     bond.SetProp("OIN_BondStereo", str(bond.GetStereo()))

        # 8. Cleanup Lone Pair Pseudos before returning
        if added_pseudo_indices and isinstance(mol, Chem.RWMol):
            PseudoAtomStrategy.remove_by_indices(mol, added_pseudo_indices)

class PseudoAtomStrategy:
    @staticmethod
    def add_lone_pair_pseudos(mol: Chem.RWMol) -> List[int]:
        """
        Adds ghost atoms to 3-coordinate P/N atoms to force RDKit chirality perception.
        """
        added_indices_and_coords = []
        conf = mol.GetConformer()
        original_num_atoms = mol.GetNumAtoms()
        
        for atom in list(mol.GetAtoms()):
            if atom.GetIdx() >= original_num_atoms: continue
            
            num_neighbors = atom.GetDegree()
            # Targeting 3rd-coordinate P/N/C that might be stereogenic
            if atom.GetAtomicNum() in (6, 7, 15) and num_neighbors == 3:
                idx = atom.GetIdx()
                pos = np.array(conf.GetAtomPosition(idx))
                
                vec_sum = np.zeros(3)
                for nbr in atom.GetNeighbors():
                    n_pos = np.array(conf.GetAtomPosition(nbr.GetIdx()))
                    vec_sum += (n_pos - pos)
                
                norm = np.linalg.norm(vec_sum)
                if norm < 0.1: continue # Nearly planar
                
                lp_vec = - (vec_sum / norm)
                lp_pos = pos + lp_vec * 1.5
                
                ghost_idx = mol.AddAtom(Chem.Atom(0)) # Wildcard *
                mol.AddBond(idx, ghost_idx, Chem.BondType.SINGLE)
                added_indices_and_coords.append((ghost_idx, lp_pos))
        
        if added_indices_and_coords:
            new_conf = Chem.Conformer(mol.GetNumAtoms())
            for i in range(original_num_atoms):
                 new_conf.SetAtomPosition(i, conf.GetAtomPosition(i))
            for g_idx, g_pos in added_indices_and_coords:
                new_conf.SetAtomPosition(g_idx, Point3D(float(g_pos[0]), float(g_pos[1]), float(g_pos[2])))
            mol.RemoveAllConformers()
            mol.AddConformer(new_conf)
            
        return [item[0] for item in added_indices_and_coords]

    @staticmethod
    def remove_by_indices(mol: Chem.RWMol, indices: List[int]) -> None:
        for idx in sorted(indices, reverse=True):
            mol.RemoveAtom(idx)

    @staticmethod
    def label_as_wildcard(mol: Chem.Mol, target_indices: List[int]) -> None:
        rw_mol = Chem.RWMol(mol)
        for idx in target_indices:
            atom = rw_mol.GetAtomWithIdx(idx)
            atom.SetAtomicNum(0)

    @staticmethod
    def strip_pseudo_atoms(mol: Chem.Mol) -> Chem.Mol:
        if not mol: return mol
        rw_mol = Chem.RWMol(mol)
        indices_to_remove = [a.GetIdx() for a in rw_mol.GetAtoms() if a.GetAtomicNum() == 0]
        for idx in sorted(indices_to_remove, reverse=True):
            rw_mol.RemoveAtom(idx)
        return rw_mol.GetMol()

class ChiralityRecoveryUtility:
    """
    Handles restoration of chiral tags and bond stereo from OIN properties.
    """
    @staticmethod
    def recover(mol: Chem.Mol):
        for atom in mol.GetAtoms():
            if atom.HasProp("OIN_CIP"):
                cip = atom.GetProp("OIN_CIP")
                if cip == "R":
                    atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                elif cip == "S":
                    atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    
    @staticmethod
    def recover_bond_stereo(mol: Chem.Mol, bond_stereo_data: List[Tuple[int, int, str]]):
        for u, v, stereo_str in bond_stereo_data:
            bond = mol.GetBondBetweenAtoms(u, v)
            if bond:
                for member in Chem.BondStereo.values.values():
                    if str(member) == stereo_str:
                        bond.SetStereo(member)
                        break
                # Handle our custom biaryl tags if RDKit doesn't recognize them
                if stereo_str.startswith("STEREOATROP"):
                    bond.SetProp("OIN_BondStereo_Recovered", stereo_str)
