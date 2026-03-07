import sys
import ast
import pprint
import argparse
from architector.complex_construction import build_complex

def main():
    parser = argparse.ArgumentParser(description='Run a single Architector job from an input dictionary file.')
    parser.add_argument('input_file', help='Path to the file containing the input dictionary (Python syntax)')
    parser.add_argument('-o', '--output', help='Base name for output XYZ files (e.g. "result"). If not provided, defaults to "output".')
    
    args = parser.parse_args()
    
    base_output_name = args.output if args.output else "output"
    
    try:
        with open(args.input_file, 'r') as f:
            content = f.read()
            
        # Parse the python dictionary safely
        input_dict = ast.literal_eval(content)
        
        print(f"Loaded input dictionary from {args.input_file}")
        
        print("\nRunning Architector build_complex...")
        output = build_complex(inputDict=input_dict)
        
        print("\nJob Completed Successfully!")
        print("Output keys:", output.keys())
        
        # Helper to find complex_molecules recursively
        def find_molecules(data):
            found = []
            if isinstance(data, dict):
                if 'complex_molecules' in data:
                    # Found the canonical list of results. Return it.
                    return data['complex_molecules']
                
                if 'ase_atoms' in data:
                    # This dictionary IS a molecule structure.
                    return [data]
                
                # Recurse into values
                for value in data.values():
                    found.extend(find_molecules(value))
            elif isinstance(data, list):
                for item in data:
                    found.extend(find_molecules(item))
            return found

        molecules = find_molecules(output)

        # Save generated structures
        if molecules:
            count = len(molecules)
            print(f"Found {count} complex molecules (recursively).")
            
            # Import ase.io inside to avoid early import errors if not needed (though needed for verify)
            try:
                from ase.io import write
            except ImportError:
                print("Warning: ASE not found. Cannot save XYZ files.")
                sys.exit(0)

            for i, mol_struct in enumerate(molecules):
                # Determine filename
                if count == 1:
                    filename = f"{base_output_name}.xyz"
                else:
                    filename = f"{base_output_name}_{i+1}.xyz"
                
                print(f"Saving molecule {i+1} to {filename}...")
                
                # Check structure type and save
                if isinstance(mol_struct, dict) and 'ase_atoms' in mol_struct:
                    write(filename, mol_struct['ase_atoms'], format='xyz')
                elif hasattr(mol_struct, 'write_xyz'): # Architector Molecule object sometimes has this
                    mol_struct.write_xyz(filename)
                elif hasattr(mol_struct, 'write_file'):
                    mol_struct.write_file(filename)
                else:
                    # Fallback or error
                    print(f"  Warning: Unknown structure format for molecule {i+1}. Skipping save.")
                    pass
        else:
            print("No 'complex_molecules' found in output (checked recursively).")
            # Print structure of top level for debugging
            print("Full Output Structure:")
            pprint.pprint(output)
        
    except FileNotFoundError:
        print(f"Error: File not found - {args.input_file}")
        sys.exit(1)
    except SyntaxError as e:
        print(f"Error: Input file is not a valid Python dictionary syntax.")
        print(e)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        # Print full traceback for debugging?
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
