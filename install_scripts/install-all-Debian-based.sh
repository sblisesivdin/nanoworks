#!/bin/bash

# Define environment and paths
ENV_NAME=".venv_nw"
INSTALL_DIR="$HOME/$ENV_NAME"
USERNAME=$(whoami)

echo "Starting installation script for GPAW and related tools..."

# Update and upgrade system packages
echo "Updating and upgrading system packages..."
sudo apt update && sudo apt upgrade -y

# Install required system packages
echo "Installing required system packages..."
sudo apt install -y python3-venv python3-pip unzip python-is-python3 \
                    python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev \
                    libfftw3-dev libkim-api-dev openkim-models libkim-api2 pkg-config \
                    task-spooler

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

# Install Nanoworks
echo "Installing Nanoworks [DFT, MD and ML]..."

# Add Nanoworks directory to PATH in .bashrc
pip install "nanoworks[all]"

echo "Creating examples folder..."
nanoworks --install-examples
echo "Examples folder is installed to ~/.nanoworks/examples ..."
# Final message
echo "Installation complete!"

