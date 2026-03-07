try:
    from openbabel import openbabel
    import sys
except ImportError:
    print("OpenBabel not found via 'openbabel'")
    try:
        import pybel
        print("Found pybel")
        ob = pybel.ob
    except ImportError:
        print("No OpenBabel found.")
        # Try finding via architector?
        import sys
        sys.exit(0)

def check_ob_indices():
    smi = "[cH]1[cH][cH][cH][cH]1"
    
    # Use OBMol
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "xyz") # or smi to smi
    
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, smi)
    
    print(f"SMILES: {smi}")
    print(f"Num Atoms: {mol.NumAtoms()}")
    
    # Add Hydrogens? Architector might add them.
    mol.AddHydrogens()
    print(f"Num Atoms (after AddHs): {mol.NumAtoms()}")
    
    for i in range(1, mol.NumAtoms() + 1): # OB is 1-based usually? Wrapper might be different.
        atom = mol.GetAtom(i)
        print(f"Atom {i-1}: {atom.GetType()} (AtomicNum: {atom.GetAtomicNum()})")
        
if __name__ == "__main__":
    check_ob_indices()
