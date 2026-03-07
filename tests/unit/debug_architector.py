from architector.io_molecule import Molecule
import sys

def check_architector_indices():
    smi = "[cH]1[cH][cH][cH][cH]1"
    
    # Architector Molecule from SMILES
    mol = Molecule(smi) # Assuming constructor takes SMILES or something similar?
    # Actually Architector usually takes inputDict.
    # But usually it uses OpenBabel internally via pybel or similar.
    # Let's try to find how it parses ligands.
    
    # If I can't instantiate Molecule directly easily, I'll simulate what build_complex does.
    # It takes 'smiles'.
    # It likely converts it to an object.
    
    print(f"SMILES: {smi}")
    try:
        # Check if Molecule has a from_smiles or acts on SMILES
        # Checking atoms
        print("Atoms in Architector Molecule:")
        for i, atom in enumerate(mol.atoms): # Assuming .atoms list
             print(f"Atom {i}: {atom.element} (Type: {atom.type})")
    except Exception as e:
        print(f"Error inspecting molecule: {e}")
        # Try inspecting mol.mol (if it wraps OB/RDKit)
        if hasattr(mol, 'mol'):
            print(f"Wrapped Mol type: {type(mol.mol)}")

if __name__ == "__main__":
    check_architector_indices()
