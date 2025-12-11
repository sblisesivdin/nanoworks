#!/bin/bash

# Define environment and paths
ENV_NAME="gpaw25_2env"
INSTALL_DIR="$HOME/$ENV_NAME"
GPAW_SETUP_DIR="$HOME/.gpaw/gpaw-setups"
USERNAME=$(whoami)

echo "Starting installation script for GPAW and related tools..."

# Update and upgrade system packages
echo "Updating and upgrading system packages..."
sudo apt update && sudo apt upgrade -y

# Install required system packages
echo "Installing required system packages..."
sudo apt install -y python3-tk python3-venv python3-pip unzip python-is-python3 \
                    python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev \
                    libfftw3-dev libkim-api-dev openkim-models libkim-api2 pkg-config \
                    task-spooler

# Create and activate the Python virtual environment
echo "Creating Python virtual environment..."
python3 -m venv "$INSTALL_DIR"
source "$INSTALL_DIR/bin/activate"

# Install ASE and other required Python packages in the virtual environment
echo "Installing ASE and required Python libraries..."
pip3 install --upgrade ase
pip3 install --upgrade setuptools_scm spglib docutils elastic requests phonopy asap3 kimpy
pip3 install --upgrade pyrapl pymongo pandas

# Set up GPAW configurations
echo "Setting up GPAW configurations..."
mkdir -p ~/.gpaw
cat > ~/.gpaw/siteconfig.py <<EOL
fftw = True
scalapack = True
libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
EOL

# Install GPAW and set up datasets
echo "Installing GPAW..."
pip3 install --upgrade gpaw

echo "Setting up PAW datasets..."
mkdir -p "$GPAW_SETUP_DIR"
gpaw install-data "$GPAW_SETUP_DIR/"

echo "Setting up Pint..."

# Add Pint directory to PATH in .bashrc
echo "Adding Pint to PATH in .bashrc..."
echo "export PATH=$HOME/pint-main:\$PATH" >> ~/.bashrc
source ~/.bashrc

# Final message
echo "Installation complete!"

