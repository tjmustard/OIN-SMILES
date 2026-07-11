import copy
import logging
import shutil

logger = logging.getLogger(__name__)

# Warn only once per process when a g-xTB optimization is requested but the
# 'xtb' binary is absent, instead of printing the fallback per conformer.
_XTB_FALLBACK_WARNED = False


class ASEOptimizer:
    """A unified wrapper for ASE-compatible calculators (xTB, MLIPs) to perform.

    geometry optimization on a given molecule object.
    """

    def __init__(self, method="xtb", fmax=0.05, max_steps=200, timeout=None):
        """Initialize the Ase optimizer."""
        self.method = method.lower()
        self.fmax = fmax
        self.max_steps = max_steps
        self.timeout = timeout

        if self.method in ("xtb", "g-xtb"):
            # We use a subprocess wrapper for g-xTB now, so no Python package imports are needed
            # here.
            pass
        elif self.method in ("mace-omol25", "mace-omol-0-extra-large-1024"):
            try:
                from mace.calculators import MACECalculator

                self._calc_cls = MACECalculator
            except ImportError:
                raise ImportError(
                    "mace-torch is not installed. Please install it via 'uv add mace-torch'."
                )
        else:
            raise ValueError(f"Optimizer method {method} not supported yet.")

    def optimize(self, mol):
        """Takes a chem.Molecule, converts to ASE Atoms, optimizes, and updates the coordinates.

        Returns a tuple of (success, energy_in_eV, optimized_mol).
        """
        # g-xTB runs as an external subprocess. If the binary is missing, return
        # the unoptimized geometry immediately -- before any deepcopy or ASE import --
        # and warn exactly once per process rather than once per conformer. This is
        # a deliberate, non-fatal degradation (contrast MACE, which raises).
        if self.method in ("xtb", "g-xtb") and shutil.which("xtb") is None:
            global _XTB_FALLBACK_WARNED
            if not _XTB_FALLBACK_WARNED:
                logger.warning(
                    "g-xTB optimizer requested but the 'xtb' binary was not found in "
                    "PATH; returning the force-field geometry with no semi-empirical "
                    "refinement. Install g-xTB to enable it. (Shown once per process.)"
                )
                _XTB_FALLBACK_WARNED = True
            return False, 0.0, mol

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
        charge = opt_mol.chg if hasattr(opt_mol, "chg") and opt_mol.chg is not None else 0
        mult = (
            opt_mol.multiplicity
            if hasattr(opt_mol, "multiplicity") and opt_mol.multiplicity is not None
            else 1
        )
        uhf = mult - 1  # number of unpaired electrons

        # Apply charge and spin using ASE's standard initializers
        charges = [0] * len(atoms)
        if charge != 0 and len(charges) > 0:
            charges[0] = charge
        atoms.set_initial_charges(charges)

        magmoms = [0] * len(atoms)
        if uhf > 0 and len(magmoms) > 0:
            magmoms[0] = uhf
        atoms.set_initial_magnetic_moments(magmoms)

        if self.method in ("xtb", "g-xtb"):
            import os
            import subprocess
            import tempfile

            from ase.io import read, write

            # The missing-binary case is handled up front in optimize(); reaching
            # here means shutil.which("xtb") succeeded.
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    input_xyz = os.path.join(tmpdir, "struc.xyz")
                    # Write the ASE atoms to XYZ
                    write(input_xyz, atoms)

                    # Create charge and uhf files if necessary
                    if charge != 0:
                        with open(os.path.join(tmpdir, ".CHRG"), "w") as f:
                            f.write(str(charge) + "\n")
                    if uhf > 0:
                        with open(os.path.join(tmpdir, ".UHF"), "w") as f:
                            f.write(str(uhf) + "\n")

                    # Run g-xTB optimization
                    # The binary supports --gxtb --opt
                    cmd = ["xtb", "struc.xyz", "--gxtb", "--opt"]
                    result = subprocess.run(
                        cmd,
                        cwd=tmpdir,
                        capture_output=True,
                        text=True,
                        check=False,
                        timeout=self.timeout,
                    )

                    if result.returncode != 0:
                        print(f"Warning: xTB failed with return code {result.returncode}.")
                        print(f"Stdout:\n{result.stdout}\nStderr:\n{result.stderr}")
                        return False, 0.0, mol

                    # Parse the optimized XYZ file
                    opt_xyz = os.path.join(tmpdir, "xtbopt.xyz")
                    if not os.path.exists(opt_xyz):
                        print(f"Warning: xTB finished but '{opt_xyz}' was not generated.")
                        return False, 0.0, mol

                    opt_atoms = read(opt_xyz)
                    atoms.set_positions(opt_atoms.get_positions())

                    # Try to parse the final energy from xtbopt.log or stdout
                    energy = 0.0
                    for line in reversed(result.stdout.splitlines()):
                        if "TOTAL ENERGY" in line:
                            parts = line.split()
                            try:
                                # xtb outputs energy in Eh (Hartrees)
                                energy_eh = float(parts[-3])
                                # Convert Hartree to eV
                                energy = energy_eh * 27.211386245988
                                break
                            except (ValueError, IndexError):
                                pass

            except subprocess.TimeoutExpired:
                print(f"Warning: xTB timed out after {self.timeout} seconds. Falling back to FF.")
                return False, 0.0, mol
            except Exception as e:
                print(f"Warning: Optimizer xTB wrapper failed. Falling back to FF. Details: {e}")
                return False, 0.0, mol

        elif self.method in ("mace-omol25", "mace-omol-0-extra-large-1024"):
            import os

            import torch

            try:
                from dotenv import load_dotenv

                load_dotenv()
            except ImportError:
                pass

            env_var = (
                "MACE_OMOL25_MODEL_PATH"
                if self.method == "mace-omol25"
                else "MACE_OMOL_0_EXTRA_LARGE_MODEL_PATH"
            )
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
                        # Test if the current GPU architecture is actually supported by this PyTorch
                        # build
                        # We must use an operation like matmul to force an actual compute kernel
                        # launch
                        _ = torch.matmul(torch.ones(10, 10).cuda(), torch.ones(10, 10).cuda())

                    import warnings

                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore",
                            message=".*TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD.*",
                            category=UserWarning,
                        )
                        calc = self._calc_cls(
                            model_paths=model_path, device=device, default_dtype="float64"
                        )

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

            # MACE attaches an ASE calculator; read the relaxed energy from it.
            # (The xtb path uses a subprocess and already parsed ``energy`` from
            # xtb's stdout above -- calling get_potential_energy() there would
            # raise "Atoms object has no calculator".)
            energy = atoms.get_potential_energy()

        new_positions = atoms.get_positions()

        for i, atom in enumerate(opt_mol.atom_list):
            atom.set_coordinate(new_positions[i].tolist())

        opt_mol.energy = energy

        return True, energy, opt_mol
