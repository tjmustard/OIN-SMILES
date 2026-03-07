#!/usr/bin/env python3
import os
import sys
import shutil
import datetime

def archive_active_specs(feature_name: str):
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    active_dir = os.path.join(base_dir, 'spec', 'active')
    archive_dir = os.path.join(base_dir, 'spec', 'archive')

    if not os.path.exists(active_dir):
        print(f"ERROR: Active directory not found at {active_dir}")
        sys.exit(1)
        
    if not os.path.exists(archive_dir):
        os.makedirs(archive_dir)

    files_to_move = [f for f in os.listdir(active_dir) if f != '.gitkeep' and os.path.isfile(os.path.join(active_dir, f))]
    
    if not files_to_move:
        print("INFO: No active specifications found to archive.")
        sys.exit(0)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    sanitized_name = "".join([c if c.isalnum() else "_" for c in feature_name])
    target_folder = os.path.join(archive_dir, f"{timestamp}_{sanitized_name}")
    
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    for filename in files_to_move:
        src = os.path.join(active_dir, filename)
        dst = os.path.join(target_folder, filename)
        shutil.move(src, dst)
        print(f"Archived: {filename} -> {target_folder}/")
        
    # Strictly enforce the .gitkeep initialized state
    gitkeep_path = os.path.join(active_dir, '.gitkeep')
    if not os.path.exists(gitkeep_path):
        with open(gitkeep_path, 'w') as f:
            pass

    print("=" * 50)
    print("ARCHIVAL SCRIPT: SUCCESS")
    print("=" * 50)
    print(f"Active directory flushed.")
    print(f"Artifacts permanently stored at absolute path:")
    print(f"{target_folder}/")
    print("=" * 50)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python archive_specs.py <Feature_Name>")
        sys.exit(1)
    archive_active_specs(sys.argv[1])