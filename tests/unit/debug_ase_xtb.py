from ase.calculators.xtb import XTB
import sys

print("Checking ASE XTB Calculator...")
try:
    # Mimic architector call
    calc = XTB(method="GFN2-xTB", solvent="water", 
               accuracy=1.0, electronic_temperature=300, max_iterations=250)
    print("[OK] ASE XTB accepted Architector arguments.")
except TypeError as e:
    print(f"[FAIL] ASE XTB rejected arguments: {e}")
    # Inspect what it accepts?
    import inspect
    print(inspect.signature(XTB.__init__))

except Exception as e:
    print(f"[ERROR] {e}")
