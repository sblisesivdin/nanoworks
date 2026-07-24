import sys
import argparse
from pathlib import Path
import os
import shutil
import importlib.resources as pkg_resources
import nanoworks

def deploy_examples():
    # Find the user's home directory
    home_dir = os.path.expanduser("~")
    target_dir = os.path.join(home_dir, ".nanoworks", "examples")
    
    # Check if the directory already exists to prevent overwriting
    if os.path.exists(target_dir):
        print(f"Warning: examples already exist in '{target_dir}'.")
        return

    print(f"Copying Nanoworks examples to '{target_dir}'...")
    
    try:
        # Locate the 'examples' directory within the installed package
        source_dir = pkg_resources.files(nanoworks).joinpath("examples")
        
        # Copy the directory tree to the target location
        shutil.copytree(source_dir, target_dir)
        print("Success! You can find the examples in the ~/.nanoworks/examples directory.")
    except Exception as e:
        print(f"An error occurred while copying the examples: {e}")

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
    parser.add_argument('-v', '--version', action='store_true', help='Show version and detailed library information')
    parser.add_argument('--install-examples', action='store_true', help='Copy example files to ~/.nanoworks/Examples')
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    # Parse the arguments. Using parse_known_args to avoid exiting on unknown arguments 
    args, unknown = parser.parse_known_args()

    # If the user passes the --install-examples flag, run the function and exit
    if args.install_examples:
        deploy_examples()
        sys.exit(0)
    
    if args.version:
        import ase
        import gpaw
        import phonopy
        try:
            import asap3
        except ImportError:
            print("----------------------------------")
            print("Welcome to Nanoworks!")
            print(f"Version: {nanoworks.__version__}")
            print("----------------------------------")
            print("Libraries used:") 
            print(f"ASE: {ase.__version__}, GPAW: {gpaw.__version__}, Phonopy: {phonopy.__version__}")
            print("----------------------------------")
            folders = ["optimizations", "examples"]
            for folder in folders:
                path = find_package_folder(folder)
                if path:
                    print(f"{folder.capitalize()} folder: {path}")
                else:
                    print(f"Could not locate {folder} folder. (It may not be included in the installation)")
            print("----------------------------------")
            print("If you do not have examples, run nanoworks --install-examples")
            print("and then continue with each example. Every example has its own README.md")
            sys.exit(1)
        
        print("----------------------------------")
        print("Welcome to Nanoworks!")
        print(f"Version: {nanoworks.__version__}")
        print("----------------------------------")
        print("Libraries used:") 
        print(f"ASE: {ase.__version__}, GPAW: {gpaw.__version__}, Phonopy: {phonopy.__version__}, ASAP3: {asap3.__version__}")
        print("----------------------------------")
        folders = ["optimizations", "examples"]
        for folder in folders:
            path = find_package_folder(folder)
            if path:
                print(f"{folder.capitalize()} folder: {path}")
            else:
                print(f"Could not locate {folder} folder. (It may not be included in the installation)")
        print("----------------------------------")
            print("If you do not have examples, run nanoworks --install-examples")
            print("and then continue with each example. Every example has its own README.md")
    

if __name__ == "__main__":
    main()
