from oinsmiles.generation.architector_adapter import ArchitectorAdapter
from oinsmiles.generation.oin_parser import OINParser
from rdkit import Chem

def debug_ferrocene():
    oin = "[Fe_LIN].[cH]{0>}1[cH]{0}[cH]{0}[cH]{0}[cH]{0}1.[cH]{1>}1[cH]{1}[cH]{1}[cH]{1}[cH]{1}1"
    parser = OINParser()
    parsed = parser.parse(oin)
    
    adapter = ArchitectorAdapter()
    res = adapter.convert(parsed)
    
    lig = res['ligands'][0]
    smi = lig['smiles']
    atoms = lig['coordinating_atoms']
    
    print(f"SMILES: {smi}")
    print(f"Binds: {atoms}")
    
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if not mol:
        print("RDKit failed to parse SMILES even with sanitize=False")
        return

    print(f"Num Atoms: {mol.GetNumAtoms()}")
    for idx, atom in enumerate(mol.GetAtoms()):
        print(f"Atom {idx}: {atom.GetSymbol()} (Map: {atom.GetAtomMapNum()})")
        
if __name__ == "__main__":
    debug_ferrocene()
