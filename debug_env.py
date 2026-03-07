import sys
import os

print(f"Python Executable: {sys.executable}")
print("System Path:")
for p in sys.path:
    print(f"  {p}")

print("\nAttempting to import architector...")
try:
    import architector
    print(f"Architector found at: {os.path.dirname(architector.__file__)}")
    print(f"Architector Version: {getattr(architector, '__version__', 'Unknown')}")
    
    print("\nAttempting submodules...")
    try:
        from architector import complex_construction
        print("SUCCESS: import architector.complex_construction")
    except ImportError as e:
        print(f"FAILURE: import architector.complex_construction: {e}")
        
    try:
        from architector import io_molecule
        print("SUCCESS: import architector.io_molecule")
    except ImportError as e:
        print(f"FAILURE: import architector.io_molecule: {e}")

except ImportError as e:
    print(f"FAILURE: Could not import architector: {e}")
