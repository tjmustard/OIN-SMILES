import sys
import traceback
import os

print(f"PYTHONPATH: {os.environ.get('PYTHONPATH')}")
print(f"Executable: {sys.executable}")

print("\nAttempting to import xtb package...")
try:
    import xtb
    print(f"[OK] xtb imported from: {xtb.__file__}")
except Exception:
    traceback.print_exc()

print("\nAttempting to import xtb.ase...")
try:
    import xtb.ase
    print(f"[OK] xtb.ase imported from: {xtb.ase.__file__}")
except Exception:
    traceback.print_exc()

print("\nAttempting to import xtb.ase.calculator...")
try:
    import xtb.ase.calculator
    print(f"[OK] xtb.ase.calculator imported from: {xtb.ase.calculator.__file__}")
except Exception:
    traceback.print_exc()

print("\nAttempting to import XTB class...")
try:
    from xtb.ase.calculator import XTB
    print(f"[OK] XTB class imported: {XTB}")
    calc = XTB(method="GFN2-xTB")
    print(f"[OK] Instantiated XTB: {calc}")
except Exception:
    traceback.print_exc()
