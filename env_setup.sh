#!/bin/bash

# FIDO Environment Setup Script
# This script installs dependencies for the "fido" conda environment

set -e  # Exit on any error

# Set installation directory
INSTALL_DIR="$(pwd)/"
BIN_DIR="$INSTALL_DIR/bin"

echo "Setting up FIDO conda environment..."
echo "Installation directory: $INSTALL_DIR"

# Create installation directories
mkdir -p "$INSTALL_DIR"
mkdir -p "$BIN_DIR"

# Create conda environment
conda create -n fido python=3.10 -y

echo "Environment 'fido' created. Activating it..."
# Note: conda activate in scripts requires special handling
eval "$(conda shell.bash hook)"
conda activate fido

echo "Installing HMMER..."
# HMMER installation
cd "$INSTALL_DIR"
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz
cd hmmer-3.4
./configure --prefix="$INSTALL_DIR/hmmer"
make
make check                        # Run automated tests
make install                      # Install HMMER programs, man pages
(cd easel; make install)          # Install Easel tools
cd "$INSTALL_DIR"

echo "Installing Clustal Omega..."
# Clustal Omega installation
cd "$INSTALL_DIR"
wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64 -O clustalo
chmod u+x clustalo
mv clustalo "$INSTALL_DIR/"

# Add installation directories to PATH for this session
export PATH="$BIN_DIR:$INSTALL_DIR/hmmer/bin:$PATH"

echo "Installing conda packages..."
# Bioinformatics tools via conda
conda install -c bioconda blast -y
conda install -c conda-forge -c bioconda mmseqs2 -y

# Python packages
conda install pandas -y
#conda install tensorflow -y
#conda install scikit-learn matplotlib -y

echo "Cleaning up..."
# Clean up downloaded files
cd "$INSTALL_DIR"
rm -f hmmer.tar.gz
rm -rf hmmer-3.4

echo "Setup complete!"
echo ""
echo "Installation directory: $INSTALL_DIR"
echo "Programs installed:"
echo "  - HMMER: $INSTALL_DIR/hmmer/"
echo "  - Clustal Omega: $BIN_DIR/clustalo"
echo ""
echo "IMPORTANT: To use this environment, run these commands:"
echo "  conda activate fido"
echo ""
echo "To make the programs permanently available, add this to your ~/.bashrc:"
echo "  export PATH=\"$BIN_DIR:$INSTALL_DIR/hmmer/bin:\$PATH\""
echo ""
echo "Then start with the 'start here for data preprocessing.ipynb' notebook"