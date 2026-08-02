#!/bin/bash

# Define environment and paths
ENV_NAME=".venv_nw"
INSTALL_DIR="$HOME/$ENV_NAME"
USERNAME=$(whoami)

echo "Starting installation script for GPAW and related tools..."
echo ""

# ---------------------------------------------------------
# INTERACTIVE INSTALLATION MENU
# ---------------------------------------------------------
echo "================================================="
echo " Select the Nanoworks installation mode: "
echo "================================================="
echo "  1) DFT only        (pip install nanoworks)"
echo "  2) DFT + MD        (pip install nanoworks[md])"
echo "  3) DFT + ML        (pip install nanoworks[ml])"
echo "  4) DFT + MD + ML   (pip install nanoworks[all]) [Default]"
echo "================================================="

# Read input directly from the terminal to prevent curl piping issues
read -r -p "Enter your choice [1-4, default=4]: " choice < /dev/tty || choice="4"

case "$choice" in
    1)
        NW_PACKAGE="nanoworks"
        MODE_NAME="DFT only"
        ;;
    2)
        NW_PACKAGE="nanoworks[md]"
        MODE_NAME="DFT + MD"
        ;;
    3)
        NW_PACKAGE="nanoworks[ml]"
        MODE_NAME="DFT + ML"
        ;;
    *)
        NW_PACKAGE="nanoworks[all]"
        MODE_NAME="DFT + MD + ML"
        ;;
esac

echo ""
echo "-> You selected: $MODE_NAME ($NW_PACKAGE)"
echo "-> Proceeding with system setup..."
echo ""
# ---------------------------------------------------------

# Update and upgrade system packages
echo "Updating and upgrading system packages..."
sudo apt update && sudo apt upgrade -y

# Install required system packages
echo "Installing required system packages..."
sudo apt install -y python3-venv python3-pip unzip python-is-python3 \
                    python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev \
                    libfftw3-dev libkim-api-dev openkim-models libkim-api2 pkg-config \
                    task-spooler build-essential

# Create and activate the Python virtual environment
echo "Creating Python virtual environment..."
python3 -m venv "$INSTALL_DIR"
source "$INSTALL_DIR/bin/activate"

# Set up GPAW configurations
echo "Setting up GPAW configurations..."
mkdir -p ~/.gpaw
cat > ~/.gpaw/siteconfig.py <<EOL
fftw = True
scalapack = True
libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
EOL

# Install Nanoworks using the user's selected package
echo "Installing $NW_PACKAGE..."
pip install "$NW_PACKAGE" --no-cache-dir

echo "Creating examples folder..."
nanoworks --install-examples
echo "Examples folder is installed to ~/.nanoworks/examples ..."

# Final message
echo "Installation complete!"
