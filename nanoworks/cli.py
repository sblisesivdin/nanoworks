import sys
import argparse
from pathlib import Path
import nanoworks

def find_package_folder(folder_name):
    """
    Attempts to locate a specific folder associated with the nanoworks package.
    Checks:
    1. Inside the package directory (e.g., site-packages/nanoworks/folder)
    2. Sibling to the package directory (e.g., repo-root/folder)
    3. sys.prefix/share/nanoworks/folder (standard data location)
    """
    # nanoworks.__file__ points to .../nanoworks/__init__.py
    package_dir = Path(nanoworks.__file__).parent
    
    # 1. Check inside package (if installed as package data)
    candidate = package_dir / folder_name
    if candidate.exists() and candidate.is_dir():
        return candidate.resolve()
        
    # 2. Check sibling (development/editable mode where folders are at repo root)
    candidate = package_dir.parent / folder_name
    if candidate.exists() and candidate.is_dir():
        return candidate.resolve()
        
    # 3. Check sys.prefix/share (system install)
    candidate = Path(sys.prefix) / "share" / "nanoworks" / folder_name
    if candidate.exists() and candidate.is_dir():
        return candidate.resolve()

    return None

def main():
    parser = argparse.ArgumentParser(prog='nanoworks', description='Nanoworks CLI tool')
    parser.add_argument('-v', '--version', action='version', version=f'nanoworks {nanoworks.__version__}')
    # Use parse_known_args to avoid exiting on unknown arguments if this is used as a wrapper
    parser.parse_args()

    print("Welcome to Nanoworks!")
    print(f"Version: {nanoworks.__version__}")
    
    folders = ["optimizations", "examples"]
    for folder in folders:
        path = find_package_folder(folder)
        if path:
            print(f"{folder.capitalize()} folder: {path}")
        else:
            print(f"Could not locate {folder} folder. (It may not be included in the installation)")

if __name__ == "__main__":
    main()