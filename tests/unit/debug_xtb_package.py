import xtb
import sys

print(f"XTB Package: {xtb}")
print(f"Version: {getattr(xtb, '__version__', 'unknown')}")
print(f"File: {getattr(xtb, '__file__', 'unknown')}")
print("Dir(xtb):")
print(dir(xtb))

try:
    from xtb import interface
    print("\nxtb.interface imported.")
except ImportError as e:
    print(f"\nImporting xtb.interface failed: {e}")

try:
    from xtb.libxtb import VER_Major
    print(f"\nxtb.libxtb available. Major Ver: {VER_Major}")
except Exception as e:
    print(f"\nxtb.libxtb check failed: {e}")
