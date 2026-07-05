import copy
import numpy as np

class ASEOptimizer:
    """
    A unified wrapper for ASE-compatible calculators (xTB, MLIPs) to perform
    geometry optimization on a given molecule object.
    """
    def __init__(self, method="xtb", fmax=0.05, max_steps=200):
        self.method = method.lower()
        self.fmax = fmax
        self.max_steps = max_steps
        
        if self.method == "xtb":
            try:
                from xtb.ase.calculator import XTB
                self._calc_cls = XTB
            except ImportError:
                raise ImportError("xtb-python is not installed. Please install it to use xTB.")
        else:
            raise ValueError(f"Optimizer method {method} not supported yet.")
            
    def optimize(self, mol):
        """
        Takes a chem.Molecule, converts to ASE Atoms, optimizes, and updates the coordinates.
        Returns a tuple of (success, energy_in_eV, optimized_mol)
        """
        try:
            from ase import Atoms
            from ase.optimize import LBFGS
        except ImportError:
            raise ImportError("ase is not installed. Please install it.")
            
        opt_mol = copy.deepcopy(mol)
        
        # Build atoms
        numbers = [atom.get_atomic_number() for atom in opt_mol.atom_list]
        positions = [atom.get_coordinate() for atom in opt_mol.atom_list]
        atoms = Atoms(numbers=numbers, positions=positions)
        
        # Parse charge and multiplicity from the MetalloGen chem.Molecule
        charge = opt_mol.chg if hasattr(opt_mol, 'chg') and opt_mol.chg is not None else 0
        mult = opt_mol.multiplicity if hasattr(opt_mol, 'multiplicity') and opt_mol.multiplicity is not None else 1
        uhf = mult - 1 # number of unpaired electrons
        
        # Apply charge and spin using ASE's standard initializers
        charges = [0] * len(atoms)
        if charge != 0 and len(charges) > 0:
            charges[0] = charge
        atoms.set_initial_charges(charges)
        
        magmoms = [0] * len(atoms)
        if uhf > 0 and len(magmoms) > 0:
            magmoms[0] = uhf
        atoms.set_initial_magnetic_moments(magmoms)
        
        if self.method == "xtb":
            calc = self._calc_cls(method="GFN2-xTB")
            atoms.calc = calc
            
        try:
            opt = LBFGS(atoms, logfile=None)
            opt.run(fmax=self.fmax, steps=self.max_steps)
            energy = atoms.get_potential_energy()
            new_positions = atoms.get_positions()
            
            for i, atom in enumerate(opt_mol.atom_list):
                atom.set_coordinate(new_positions[i].tolist())
                
            opt_mol.energy = energy
            
            return True, energy, opt_mol
        except Exception as e:
            print(f"ASE Optimization ({self.method}) failed: {e}")
            return False, float('inf'), mol
