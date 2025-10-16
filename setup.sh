#!/bin/bash
# setup.sh - One-time setup for MetaMAG Explorer

echo "Setting up MetaMAG Explorer..."

# Make scripts executable
chmod +x install_metamag.sh
chmod +x metamag

# Create symbolic link in conda environment
CONDA_PREFIX="${CONDA_PREFIX:-$HOME/miniconda3}"
mkdir -p "$CONDA_PREFIX/bin"
ln -sf "$PWD/metamag" "$CONDA_PREFIX/bin/metamag"

echo " Setup complete!"
echo ""
echo "Next steps:"
echo "1. Run: ./install_metamag.sh"
echo "2. After installation, you can use: metamag --help"
