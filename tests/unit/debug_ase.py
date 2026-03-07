import ase
import ase.calculators
import sys

print(f"ASE Version: {ase.__version__}")
print(f"ASE Path: {ase.__file__}")

# Check for xtb in calculators
try:
    import ase.calculators.xtb
    print("Found ase.calculators.xtb")
except ImportError:
    print("Could NOT find ase.calculators.xtb")

# Check what IS in calculators
print("ase.calculators dir:", dir(ase.calculators))

# Check for generic calculators
try:
    from ase.calculators.calculator import Calculator, FileIOCalculator
    print("Found generic Calculator classes")
except ImportError:
    print("Could not find generic Calculator classes")
