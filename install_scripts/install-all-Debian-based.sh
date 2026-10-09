#!/bin/bash
# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

set -euo pipefail

# Release assets set this before invoking the installer.
NANOWORKS_VERSION="${NANOWORKS_VERSION:-}"
INSTALL_DIR="$HOME/.venv_nw"
dry_run=false
case "${1:-}" in
    --dry-run) dry_run=true ;;
    --help)
        echo "Usage: bash install-all-Debian-based.sh [--dry-run]"
        echo "NANOWORKS_ENGINE=qe|gpaw|both (default: both)"
        echo "NANOWORKS_COMPONENTS=dft|md|ml|all (default: all)"
        echo "NANOWORKS_VERSION optionally pins the Python package version."
        exit 0 ;;
    '') ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
esac
if (( $# > 1 )); then
    echo "Only one installer option is accepted." >&2
    exit 2
fi
if [[ -n "$NANOWORKS_VERSION" && ! "$NANOWORKS_VERSION" =~ ^[0-9]+(\.[0-9]+){2}([A-Za-z0-9.-]+)?$ ]]; then
    echo "Invalid NANOWORKS_VERSION: $NANOWORKS_VERSION" >&2
    exit 2
fi

echo "Starting Nanoworks installation..."
if [[ -z "${NANOWORKS_ENGINE:-}" ]]; then
    echo "DFT engine: 1) QE  2) GPAW  3) Both [default]"
    engine_choice="3"
    if [[ -t 2 ]]; then
        read -r -p "Choose engine [1-3, default=3]: " engine_choice < /dev/tty || engine_choice="3"
    fi
    case "$engine_choice" in
        1) NANOWORKS_ENGINE=qe ;;
        2) NANOWORKS_ENGINE=gpaw ;;
        3|'') NANOWORKS_ENGINE=both ;;
        *) echo "Invalid engine choice: $engine_choice" >&2; exit 2 ;;
    esac
fi
if [[ -z "${NANOWORKS_COMPONENTS:-}" ]]; then
    echo "Components: 1) DFT only  2) DFT + MD  3) DFT + ML  4) DFT + MD + ML [default]"
    component_choice="4"
    if [[ -t 2 ]]; then
        read -r -p "Choose components [1-4, default=4]: " component_choice < /dev/tty || component_choice="4"
    fi
    case "$component_choice" in
        1) NANOWORKS_COMPONENTS=dft ;;
        2) NANOWORKS_COMPONENTS=md ;;
        3) NANOWORKS_COMPONENTS=ml ;;
        4|'') NANOWORKS_COMPONENTS=all ;;
        *) echo "Invalid component choice: $component_choice" >&2; exit 2 ;;
    esac
fi

# Build the package extras independently: [all] would always pull in GPAW.
extras=()
system_packages=(python3-venv python3-pip unzip python-is-python3 task-spooler)
use_gpaw=false
use_qe=false
case "$NANOWORKS_ENGINE" in
    qe) extras+=(qe); use_qe=true ;;
    gpaw) extras+=(gpaw); use_gpaw=true ;;
    both) extras+=(gpaw qe); use_gpaw=true; use_qe=true ;;
    *) echo "Invalid NANOWORKS_ENGINE: $NANOWORKS_ENGINE" >&2; exit 2 ;;
esac
if "$use_gpaw"; then
    system_packages+=(python3-dev build-essential libopenblas-dev libxc-dev
                      libscalapack-mpi-dev libfftw3-dev pkg-config)
fi
case "$NANOWORKS_COMPONENTS" in
    dft) ;;
    md|ml|all)
        if [[ "$NANOWORKS_COMPONENTS" == md || "$NANOWORKS_COMPONENTS" == all ]]; then
            extras+=(md)
            system_packages+=(python3-dev build-essential libopenblas-dev pkg-config
                              libkim-api-dev openkim-models libkim-api2 lammps)
        fi
        if [[ "$NANOWORKS_COMPONENTS" == ml || "$NANOWORKS_COMPONENTS" == all ]]; then
            extras+=(ml)
            system_packages+=(python3-dev build-essential)
        fi ;;
    *) echo "Invalid NANOWORKS_COMPONENTS: $NANOWORKS_COMPONENTS" >&2; exit 2 ;;
esac
extra_list=$(IFS=,; echo "${extras[*]}")
NW_PACKAGE="nanoworks[${extra_list}]"
if [[ -n "$NANOWORKS_VERSION" ]]; then
    NW_PACKAGE="${NW_PACKAGE}==${NANOWORKS_VERSION}"
fi

echo "Engine: $NANOWORKS_ENGINE; components: $NANOWORKS_COMPONENTS"
echo "Python package: $NW_PACKAGE"
echo "System packages: ${system_packages[*]}"
echo "Virtual environment: $INSTALL_DIR"
echo "GPAW configuration: $use_gpaw"
echo "QE pseudopotential installation: $use_qe"
echo "QE calculations use externally installed Quantum ESPRESSO; elasticity also needs thermo_pw."
if "$dry_run"; then
    echo "Dry run complete; no installation commands executed."
    exit 0
fi

sudo apt update
sudo apt install -y "${system_packages[@]}"
python3 -m venv "$INSTALL_DIR"
source "$INSTALL_DIR/bin/activate"

if "$use_gpaw"; then
    mkdir -p "$HOME/.gpaw"
    # Preserve an existing, potentially site-specific build configuration.
    if [[ ! -e "$HOME/.gpaw/siteconfig.py" ]]; then
        cat > "$HOME/.gpaw/siteconfig.py" <<'CONFIG'
fftw = True
scalapack = True
libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
CONFIG
    fi
fi
python -m pip install "$NW_PACKAGE" --no-cache-dir
nanoworks --install-examples
if "$use_qe"; then
    nanoworks --install-qe-pseudos
fi

echo "Installation complete. Activate with: source $INSTALL_DIR/bin/activate"
echo "Examples: ~/.nanoworks/examples"
if "$use_qe"; then
    echo "QE pseudopotentials: ~/.nanoworks/pseudos/qe/pseudodojo/pbe/"
    echo "Provide QE 7.4.1 in PATH; elasticity also requires thermo_pw 2.1.0."
fi
