
from rdkit import Chem

def test_frag_construction():
    print("Test: Aromatic Bonds but IsAromatic=False on Atoms")
    
    mw = Chem.RWMol()
    # Create 6 carbons, loop
    indices = []
    for i in range(6):
        a = Chem.Atom(6) # C
        a.SetIsAromatic(False) # Explicitly False (default)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        # xyz2mol sets explicit H based on neighbors. For benzene-like radical:
        # 5 atoms have 1 H, 1 atom has 0 H.
        if i == 0:
            a.SetNumExplicitHs(0) # Radical site
        else:
            a.SetNumExplicitHs(1)
            
        indices.append(mw.AddAtom(a))
        
    # Add Aromatic Bonds
    for i in range(6):
        mw.AddBond(indices[i], indices[(i+1)%6], Chem.BondType.AROMATIC)
        
    kwmol = mw.GetMol()
    
    # Simulate OINSanitizer logic on Atom 0
    a0 = kwmol.GetAtomWithIdx(0)
    print(f"Atom 0: IsAro={a0.GetIsAromatic()}")
    
    # Calculate Deficit
    # Valence 4.
    # Bond Sum (AROMATIC=1.5). 1.5+1.5 = 3.
    # Explicit H = 0.
    # Deficit = 1.
    a0.SetNumRadicalElectrons(1)
    
    # Generate SMILES
    try:
        kwmol.UpdatePropertyCache(strict=False)
    except:
        pass
        
    s = Chem.MolToSmiles(kwmol, isomericSmiles=True, canonical=True)
    print(f"Result: {s}")
    
    if "[cH]" in s:
        print("FAIL: Got [cH]")
    elif "[c]" in s:
        print("PASS: Got [c]")
    else:
        print(f"Other: {s}")

if __name__ == "__main__":
    test_frag_construction()
