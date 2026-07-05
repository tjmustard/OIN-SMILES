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
        elif self.method in ("mace-omol25", "mace-omol-0-extra-large-1024"):
            try:
                from mace.calculators import MACECalculator
                self._calc_cls = MACECalculator
            except ImportError:
                raise ImportError("mace-torch is not installed. Please install it via 'uv add mace-torch'.")
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
            try:
                calc = self._calc_cls(method="GFN2-xTB")
                atoms.calc = calc
                opt = LBFGS(atoms, logfile=None)
                opt.run(fmax=self.fmax, steps=self.max_steps)
            except Exception as e:
                print(f"Warning: Optimizer xTB failed. Falling back to FF. Details: {e}")
                return False, 0.0, mol
                
        elif self.method in ("mace-omol25", "mace-omol-0-extra-large-1024"):
            import os
            import torch
            
            try:
                from dotenv import load_dotenv
                load_dotenv()
            except ImportError:
                pass
            
            env_var = "MACE_OMOL25_MODEL_PATH" if self.method == "mace-omol25" else "MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH"
            model_path = os.environ.get(env_var)
            
            if not model_path or not os.path.exists(model_path):
                raise ValueError(
                    f"{env_var} environment variable is not set or file does not exist. "
                    f"Please set it to the path of the {self.method} .model weights file."
                )
            
            devices_to_try = []
            if torch.cuda.is_available():
                devices_to_try.append("cuda")
            devices_to_try.append("cpu")
            
            success = False
            for device in devices_to_try:
                try:
                    # Reset atoms positions in case a previous run failed halfway through
                    atoms.set_positions(positions)
                    
                    if device == "cuda":
                        # Test if the current GPU architecture is actually supported by this PyTorch build
                        # We must use an operation like matmul to force an actual compute kernel launch
                        _ = torch.matmul(torch.ones(10, 10).cuda(), torch.ones(10, 10).cuda())
                        
                    import warnings
                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore", 
                            message=".*TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD.*", 
                            category=UserWarning
                        )
                        calc = self._calc_cls(model_paths=model_path, device=device, default_dtype="float64")
                    
                    atoms.calc = calc
                    opt = LBFGS(atoms, logfile=None)
                    opt.run(fmax=self.fmax, steps=self.max_steps)
                    success = True
                    break
                except Exception as e:
                    print(f"Warning: MACE failed on device '{device}'. Details: {e}")
                    if device == "cuda":
                        print("Falling back to CPU...")
                    else:
                        print("Falling back to FF...")
            
            if not success:
                return False, 0.0, mol
                
        energy = atoms.get_potential_energy()
        new_positions = atoms.get_positions()
        
        for i, atom in enumerate(opt_mol.atom_list):
            atom.set_coordinate(new_positions[i].tolist())
            
        opt_mol.energy = energy
        
        return True, energy, opt_mol
