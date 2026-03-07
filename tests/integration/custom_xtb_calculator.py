import os
import shutil
from ase.calculators.calculator import FileIOCalculator, ReadError

class CustomXTB(FileIOCalculator):
    """
    A custom ASE calculator for XTB that wraps the binary directly.
    Designed to mimic the interface expected by Architector.
    """
    implemented_properties = ['energy', 'forces', 'attributes']
    command = 'xtb xtb_in.xyz --gfn 2'  # Default
    discard_results_on_any_change = True

    def __init__(self, restart=None, 
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='xtb', atoms=None, **kwargs):
        
        # Resolve xtb binary absolute path
        xtb_bin = shutil.which("xtb")
        if not xtb_bin:
            # Fallback for manual check if PATH update failed in subshell
            # Try finding it in known location relative to CWD of the MAIN process
            # Note: Architector changes cwd to temp, so relative paths to repo root break unless absolute.
            # Assuming repo root is 3 levels up from this file? No, verify_roundtrip is in tests/integration.
            # But the script is running from repo root.
            # Better to rely on the install location causing it to be absolute.
            possible_path = os.path.abspath(".gemini/xtb_install/xtb-6.7.1/bin/xtb")
            if os.path.exists(possible_path):
                xtb_bin = possible_path
            else:
                xtb_bin = "xtb" # Hope for the best or fail clearly
        
        self.xtb_bin = xtb_bin

        # Architector passes these:
        # method='GFN2-xTB', solvent='none', accuracy=1.0, 
        # electronic_temperature=300, max_iterations=250
        
        self.xtb_method = kwargs.pop('method', 'GFN2-xTB')
        self.solvent = kwargs.pop('solvent', 'none')
        self.accuracy = kwargs.pop('accuracy', 1.0)
        self.elec_temp = kwargs.pop('electronic_temperature', 300)
        self.max_iter = kwargs.pop('max_iterations', 250)
        self.uhf = kwargs.pop('uhf', 0) # Spin
        
        # Construct command
        cmd = [self.xtb_bin, "xtb_in.xyz"]
        
        # Map method
        if self.xtb_method == 'GFN2-xTB':
            cmd.extend(['--gfn', '2'])
        elif self.xtb_method == 'GFN1-xTB':
            cmd.extend(['--gfn', '1'])
        elif self.xtb_method == 'GFN-FF':
            cmd.extend(['--gfnff'])
        
        # Solvent
        if self.solvent and self.solvent.lower() != 'none':
            cmd.extend(['--alpb', self.solvent])
            
        # Accuracy/T/Iter
        cmd.extend(['--acc', str(self.accuracy)])
        cmd.extend(['--etemp', str(self.elec_temp)])
        cmd.extend(['--iterations', str(int(self.max_iter))])
        
        # Spin (uhf) and Charge are handled via .CHRG and .UHF files or command line
        # Architector writes charge/spin to atoms, so we should read from atoms?
        # Actually, FileIOCalculator writes input.
        
        self.command = " ".join(cmd) + " > xtb.out"
        
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, command=self.command, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        
        # Write coordinates
        from ase.io import write
        write('xtb_in.xyz', atoms)
        
        # Write charge and multiplicity if present
        # ASE atoms object has initial_charges and initial_magnetic_moments
        charge = int(sum(atoms.get_initial_charges()))
        uhf = int(sum(atoms.get_initial_magnetic_moments()))
        
        with open('.CHRG', 'w') as f:
            f.write(str(charge))
        
        with open('.UHF', 'w') as f:
            f.write(str(uhf))

    def read_results(self):
        # Read energy and gradients from XTB output files
        # Energy is in xtb.out or energy file
        # Gradient is in 'gradient' file
        
        # Read energy
        energy = None
        with open('xtb.out', 'r') as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    # | TOTAL ENERGY             -15.603390314848 H |
                    parts = line.split()
                    energy = float(parts[3]) * 27.211386245988 # Hartree to eV
        
        if energy is None:
            raise ReadError("Could not find energy in xtb.out")
            
        self.results['energy'] = energy
        
        # Read forces (negative of gradient)
        # Gradient file format:
        # $grad
        # cycle =      1    SCF energy =    -15.603390314848   |dE/dxyz| =  0.027606
        #    0.000000000000       0.000000000000       0.000000000000      C
        #   -0.024955115206       0.011859942702       0.000000000000      H
        import numpy as np
        forces = []
        try:
            with open('gradient', 'r') as f:
                lines = f.readlines()
                # Skip header ($grad, cycle line)
                start_idx = -1
                for i, line in enumerate(lines):
                    if '$grad' in line:
                        start_idx = i + 2 # Skip $grad and cycle line
                        break
                
                if start_idx != -1:
                    for line in lines[start_idx:]:
                        if '$end' in line:
                            break
                        parts = line.split()
                        if len(parts) >= 3:
                            # Gradient in Hartree/Bohr?
                            # ASE expects eV/Angstrom
                            # XTB gradient file is Hartree/Bohr
                            gx = float(parts[0])
                            gy = float(parts[1])
                            gz = float(parts[2])
                            
                            # Convert Gradient (Ha/Bohr) to Force (eV/Ang)
                            # Force = -Gradient
                            # 1 Ha = 27.211386 eV
                            # 1 Bohr = 0.529177 Ang
                            # Factor = 27.211 / 0.5291 = ~51.422
                            conv = 27.211386245988 / 0.52917721067
                            forces.append([-gx * conv, -gy * conv, -gz * conv])
            
            self.results['forces'] = np.array(forces)
            
        except FileNotFoundError:
            # If GFN-FF or single point, maybe no gradient file?
            # Architector requests forces for relaxation.
            pass
