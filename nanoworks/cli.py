import sys
import argparse
from importlib import metadata
from pathlib import Path
import os
import shutil
import importlib.resources as pkg_resources
import nanoworks
from nanoworks.pseudos import install_qe_pseudopotentials


def format_qe_pseudo_installation(results):
    """Render a concise user-facing pseudopotential installation report."""
    scalar = results['scalar']
    family_value = str(scalar.get('family', 'pseudodojo')).lower()
    family = (
        'PseudoDojo'
        if family_value == 'pseudodojo'
        else family_value.title()
    )
    xc = str(scalar.get('xc', 'pbe')).upper()
    file_format = str(scalar.get('format', 'upf')).upper()
    accuracy = str(scalar.get('accuracy', 'standard'))

    lines = [
        'Quantum ESPRESSO pseudopotentials are ready.',
        'Library: {} {} norm-conserving {} sets'.format(
            family,
            xc,
            file_format,
        ),
        'Accuracy profile: ' + accuracy,
        'Installed sets:',
    ]

    labels = {
        'scalar': 'Scalar relativistic',
        'full': 'Fully relativistic',
    }
    for name in ('scalar', 'full'):
        result = results[name]
        status = 'already present' if result['skipped'] else 'installed'
        lines.extend([
            '  {}:'.format(labels[name]),
            '    Status: {}'.format(status),
            '    Version: {}'.format(result.get('version', 'unknown')),
            '    Files: {}'.format(result.get('count', 'unknown')),
            '    Table: {}'.format(result.get('table', 'unknown')),
            '    Directory: {}'.format(result['directory']),
            '    Manifest: {}'.format(result['manifest']),
        ])

    lines.extend([
        'Default use:',
        '  Scalar-relativistic set: standard Nanoworks QE calculations',
        '  Fully-relativistic set: available for SOC workflows',
        'Next: validate an input with dftsolve --check or '
        'dftconverge --check.',
    ])
    return '\n'.join(lines)


def _installed_version(distribution_name, optional=False):
    """Return an installed distribution version without importing it."""
    try:
        return metadata.version(distribution_name)
    except metadata.PackageNotFoundError:
        suffix = " (optional)" if optional else ""
        return "not installed" + suffix


def print_version_information():
    """Print package versions for both minimal and extended installs."""
    dependencies = (
        ("ASE", "ase", False),
        ("GPAW", "gpaw", True),
        ("Phonopy", "phonopy", True),
        ("ASAP3", "asap3", True),
    )

    print("--------------------------------------------------------------------")
    print("Welcome to Nanoworks!")
    print(f"Version: {nanoworks.__version__}")
    print("--------------------------------------------------------------------")
    print("Libraries used:")
    for label, distribution_name, optional in dependencies:
        print(
            f"{label}: "
            f"{_installed_version(distribution_name, optional=optional)}"
        )
    print("--------------------------------------------------------------------")

    examples_path = find_package_folder("examples")
    if examples_path:
        print(f"Examples folder: {examples_path}")
    else:
        print(
            "Could not locate examples folder. "
            "(It may not be included in the installation)"
        )

    print("--------------------------------------------------------------------")
    print("If you do not have examples, run nanoworks --install-examples")
    print("and then continue with each example. Every example has its own README.md")
    print("You can install pseudopotantions for QE with --install-qe-pseudos")

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
    parser.add_argument(
    '--install-qe-pseudos',
    action='store_true',
    help=(
        'Install the default Quantum ESPRESSO '
        'pseudopotential library'
    ),
)
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    # Parse the arguments. Using parse_known_args to avoid exiting on unknown arguments 
    args, unknown = parser.parse_known_args()

    # If the user passes the --install-examples flag, run the function and exit
    if args.install_examples:
        deploy_examples()
        sys.exit(0)
    
    if args.install_qe_pseudos:
        print('Installing Quantum ESPRESSO pseudopotentials...')

        results = install_qe_pseudopotentials()
        print(format_qe_pseudo_installation(results))

        sys.exit(0)
    
    if args.version:
        print_version_information()
    

if __name__ == "__main__":
    main()
