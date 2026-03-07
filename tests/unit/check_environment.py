import shutil
import sys
import os
import subprocess

def check_xtb_binary():
    print("Checking for XTB binary...")
    xtb_path = shutil.which("xtb")
    if xtb_path:
        print(f"[OK] Found XTB binary at: {xtb_path}")
        # Try version
        try:
            res = subprocess.run([xtb_path, "--version"], capture_output=True, text=True)
            print(f"    Version info: {res.stdout.splitlines()[0] if res.stdout else 'Unknown'}")
            return True
        except Exception as e:
            print(f"[WARN] Found binary but failed to run: {e}")
            return False
    else:
        print("[FAIL] XTB binary NOT found in PATH.")
        # Check standard locations?
        return False

def check_xtb_python():
    print("\nChecking for XTB Python package...")
    try:
        import xtb
        print(f"[OK] xtb package imported. Version: {getattr(xtb, '__version__', 'unknown')}")
        
        # Check if interface is usable (detects C-API failure)
        try:
            from xtb import interface
            print("    [OK] xtb.interface imported (C-API available).")
            return True
        except ImportError as e:
            print(f"    [FAIL] xtb.interface failed to import: {e}")
            print("           This usually means the C extension is missing or incompatible.")
            return False
            
    except ImportError:
        print("[FAIL] xtb python package not found.")
        return False
    except Exception as e:
        print(f"[WARN] xtb imported but error: {e}")
        return False

def main():
    print("=== OIN Engine Environment Check ===\n")
    bin_ok = check_xtb_binary()
    py_ok = check_xtb_python()
    
    if not bin_ok and not py_ok:
        print("\n[CRITICAL] XTB is missing completely. Architector cannot run GFN2-xTB.")
        print("Please install 'xtb' binary (e.g. conda install -c conda-forge xtb) or ensure it is in PATH.")
        sys.exit(1)
    
    if not bin_ok and py_ok:
        print("\n[WARN] XTB Python package found, but 'xtb' binary missing.")
        print("If Architector relies on the binary shell command, it will fail.")
        
    print("\nEnvironment check completed.")

if __name__ == "__main__":
    main()
