"""
Molassembler Adapter for 3D structure generation.
Replaces Architector with SCINE Molassembler.
"""
import os
from typing import Dict, Any, List, Optional
import numpy as np
from rdkit import Chem

from .oin_parser import ParsedOIN, OINVector

class MolassemblerTimeoutError(RuntimeError):
    """Raised when Molassembler distance geometry generation times out."""
    pass

def _molassembler_worker(args: dict) -> str:
    """
    Picklable worker function for Molassembler structure generation.
    Executed in a separate process for timeout safety.
    """
    import scine_molassembler as masm
    import scine_utilities as utils
    
    try:
        # 1. Initialize Metal
        metal_symbol = args['metal_symbol']
        try:
            metal_type = getattr(utils.ElementType, metal_symbol)
        except AttributeError:
            return f"Error: Invalid metal element type '{metal_symbol}'"
            
        mol = masm.Molecule(metal_type)
        metal_idx = 0
        
        # 2. Add Ligands
        site_to_slot = {}
        site_counter = 0
        
        for ligand in args['ligands']:
            lig_smiles = ligand['smiles']
            bond_requests = ligand['bonds'] # List of (rd_idx, slot_idx)
            
            # Use unsanitized parsing for OIN fragments
            rd_mol = Chem.MolFromSmiles(lig_smiles, sanitize=False)
            if not rd_mol:
                return f"Error: Could not parse ligand SMILES '{lig_smiles}'"
            
            rd_mol.UpdatePropertyCache(strict=False)
            rd_mol = Chem.AddHs(rd_mol)
            
            # Map rdkit_idx -> masm_idx
            rd_to_masm = {}
            coord_atoms_map = {r_idx: s_idx for r_idx, s_idx in bond_requests}
            coord_atoms = list(coord_atoms_map.keys())
            
            # BFS starting from coordination atoms.
            queue = []
            # We must be careful to add the FIRST connection to metal only once.
            first_r_idx = coord_atoms[0] if coord_atoms else None
            
            if first_r_idx is not None:
                sym = rd_mol.GetAtomWithIdx(first_r_idx).GetSymbol()
                m_type = getattr(utils.ElementType, sym)
                # Adds atom and a bond to metal. This is Site #site_counter.
                m_idx = mol.add_atom(m_type, metal_idx, masm.BondType.Single)
                rd_to_masm[first_r_idx] = m_idx
                site_to_slot[site_counter] = coord_atoms_map[first_r_idx]
                site_counter += 1
                queue.append(first_r_idx)
            
            visited = set(rd_to_masm.keys())
            head = 0
            while head < len(queue):
                curr_rd_idx = queue[head]
                head += 1
                curr_masm_idx = rd_to_masm[curr_rd_idx]
                
                for neighbor in rd_mol.GetAtomWithIdx(curr_rd_idx).GetNeighbors():
                    n_rd_idx = neighbor.GetIdx()
                    if n_rd_idx not in visited:
                        n_sym = neighbor.GetSymbol()
                        n_type = getattr(utils.ElementType, n_sym)
                        n_masm_idx = mol.add_atom(n_type, curr_masm_idx, masm.BondType.Single)
                        rd_to_masm[n_rd_idx] = n_masm_idx
                        visited.add(n_rd_idx)
                        queue.append(n_rd_idx)
                    else:
                        n_masm_idx = rd_to_masm[n_rd_idx]
                        if not mol.graph.adjacent(curr_masm_idx, n_masm_idx):
                            mol.add_bond(curr_masm_idx, n_masm_idx, masm.BondType.Single)
            
            # Add remaining bonds from multidentate ligand to metal
            for r_idx, s_idx in bond_requests:
                m_idx = rd_to_masm.get(r_idx)
                if m_idx is not None and not mol.graph.adjacent(metal_idx, m_idx):
                    # This adds a bond to an existing atom. This is Site #site_counter.
                    mol.add_bond(metal_idx, m_idx, masm.BondType.Single)
                    site_to_slot[site_counter] = s_idx
                    site_counter += 1
        
        # 3. Set Geometry
        shape_map = {
            'LIN': masm.shapes.Shape.Line,
            'SPL': masm.shapes.Shape.Square,
            'OCT': masm.shapes.Shape.Octahedron,
            'TET': masm.shapes.Shape.Tetrahedron,
            'TBP': masm.shapes.Shape.TrigonalBipyramid,
            'TPY': masm.shapes.Shape.TrigonalPyramid,
            'SPY': masm.shapes.Shape.SquarePyramid,
            'TPL': masm.shapes.Shape.EquilateralTriangle,
            'PBP': masm.shapes.Shape.PentagonalBipyramid
        }
        
        geo_code = args['geometry_code']
        shape = shape_map.get(geo_code, masm.shapes.Shape.Square)
        mol.set_shape_at_atom(metal_idx, shape)

        # 3.5. Enforce OIN Slot Mapping
        p = mol.stereopermutators.option(metal_idx)
        if p and p.num_assignments > 1:
            best_assignment = None
            shape_coords = masm.shapes.coordinates(shape)
            
            if len(site_to_slot) >= 2:
                for assignment in range(p.num_assignments):
                    mol.assign_stereopermutator(metal_idx, assignment)
                    current_p = mol.stereopermutators.option(metal_idx)
                    
                    match = True
                    # Check all possible pairs of sites to be exhaustive
                    for i in range(site_counter):
                        for j in range(i + 1, site_counter):
                            if i in site_to_slot and j in site_to_slot:
                                target_angle = np.arccos(np.clip(
                                    np.dot(shape_coords[:, site_to_slot[i]], shape_coords[:, site_to_slot[j]]), 
                                    -1.0, 1.0))
                                actual_angle = current_p.angle(i, j)
                                if not np.isclose(target_angle, actual_angle, atol=0.1):
                                    match = False
                                    break
                        if not match: break
                    
                    if match:
                        best_assignment = assignment
                        break
                
                if best_assignment is not None:
                    mol.assign_stereopermutator(metal_idx, best_assignment)
        
        # 4. Generate Conformation
        config = masm.dg.Configuration()
        res = masm.dg.generate_conformation(mol, args['seed'], config)
        
        if isinstance(res, masm.dg.Error):
            return f"Error: Molassembler DG Failure: {res.message}"
        
        if res is None:
            return "Error: Molassembler DG returned None"
            
        # 6. Extract XYZ
        positions = res
        # Molassembler DG output is in Bohr. Convert to Angstroms.
        BOHR_TO_ANGSTROM = 0.5291772109
        positions *= BOHR_TO_ANGSTROM
        
        output = []
        output.append(str(mol.graph.V))
        output.append(f"Generated by OIN-SMILES MolassemblerAdapter | Geo: {geo_code}")
        
        for i in range(mol.graph.V):
            elem = mol.graph.element_type(i)
            elem_str = str(elem).split('.')[-1]
            pos = positions[i]
            output.append(f"{elem_str} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}")
            
        return "\n".join(output)
        
    except Exception as e:
        import traceback
        return f"Error: Unexpected exception in worker:\n{traceback.format_exc()}"

class MolassemblerAdapter:
    def __init__(self, timeout: int = 60, seed: int = 42):
        self.timeout = timeout
        self.seed = seed

    def generate(self, parsed_oin: ParsedOIN) -> str:
        """
        Orchestrates 3D generation via Molassembler in a separate process.
        """
        # 1. Prepare Arguments
        metal_frag = parsed_oin.fragments[parsed_oin.metal_fragment_idx]
        metal_symbol = self._get_metal_symbol(metal_frag)
        
        ligands = []
        # Group vectors by fragment
        frag_data = {}
        for v in parsed_oin.vectors:
            if v.fragment_idx not in frag_data:
                frag_data[v.fragment_idx] = []
            frag_data[v.fragment_idx].append((v.atom_in_fragment_idx, v.slot_idx))
            
        for i, frag_smiles in enumerate(parsed_oin.fragments):
            if i == parsed_oin.metal_fragment_idx:
                continue
            
            bonds = frag_data.get(i, [])
            if not bonds:
                continue
                
            ligands.append({
                'smiles': frag_smiles,
                'bonds': bonds
            })
            
        args = {
            'metal_symbol': metal_symbol,
            'ligands': ligands,
            'geometry_code': parsed_oin.geometry,
            'seed': self.seed
        }
        
        # 2. Run in Process via ProcessPoolExecutor as mandated by MiniPRD
        import concurrent.futures
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
            future = executor.submit(_molassembler_worker, args)
            try:
                result = future.result(timeout=self.timeout)
            except concurrent.futures.TimeoutError:
                raise MolassemblerTimeoutError(f"Molassembler generation timed out after {self.timeout}s")
            except Exception as e:
                # Catching any worker-level exceptions that bubble up
                raise RuntimeError(f"Error in Molassembler worker: {str(e)}")
            
        if result.startswith("Error:"):
            raise RuntimeError(result)
            
        return result

    def _get_metal_symbol(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumAtoms() == 1:
            return mol.GetAtomWithIdx(0).GetSymbol()
        return smiles.replace("[", "").replace("]", "")
