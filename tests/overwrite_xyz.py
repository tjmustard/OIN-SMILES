import sys
import re
import os

def overwrite_with_xyz(filename):
    # 1. Read the original data into memory
    try:
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return

    # 2. Process the data
    atom_count = len(lines)
    output_lines = []

    # XYZ Header
    output_lines.append(str(atom_count))
    output_lines.append("") # Blank comment line

    for line in lines:
        parts = line.split()
        
        # Ensure line has enough columns before processing
        if len(parts) < 5:
            print(f"Skipping malformed line in {filename}: {line}")
            continue

        raw_label = parts[1] # e.g., Pt1
        x, y, z = parts[2], parts[3], parts[4]

        # Remove digits from the element label (Pt1 -> Pt)
        element = re.sub(r'\d+', '', raw_label)

        # Format: Left-align element, Right-align coords
        formatted_line = f"{element:<2} {x:>10} {y:>10} {z:>10}"
        output_lines.append(formatted_line)

    # 3. Overwrite the original file
    try:
        with open(filename, 'w') as f:
            f.write("\n".join(output_lines))
            f.write("\n")
        
        print(f"Overwrote {filename} with valid XYZ format ({atom_count} atoms).")

    except Exception as e:
        print(f"Error writing to file: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python overwrite_xyz.py <filename>")
    else:
        # Loop through all arguments so you can pass multiple files if needed
        # e.g., python overwrite_xyz.py file1.txt file2.txt
        for target_file in sys.argv[1:]:
            overwrite_with_xyz(target_file)