import sys
import traceback

print("Attempting to import xtb.interface...")
try:
    import xtb.interface
    print("[OK] xtb.interface imported successfully.")
except Exception:
    traceback.print_exc()

print("\nAttempting to import xtb.ase.calculator...")
try:
    from xtb.ase.calculator import XTB
    print("[OK] xtb.ase.calculator.XTB imported successfully.")
except Exception:
    traceback.print_exc()
