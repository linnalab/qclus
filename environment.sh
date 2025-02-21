#!/bin/bash
set -e # Exits immediately if a command exits with a non-zero status

# Step 1: Create the Conda environment
echo "Creating the Conda environment..."
conda create --name qclus python=3.11 -y

# Step 2: Activate the environment
echo "Activating the environment..."
source activate qclus

# Step 3: Install packages using Conda
echo "Installing Conda packages..."
conda install -c conda-forge scanpy python-igraph leidenalg loompy -y
conda install -c conda-forge python-annoy -y
conda install -c conda-forge jupyterlab -y

# Step 4: Install additional Python packages using pip
echo "Installing Python packages using pip..."
pip install .

# Joke sequence
echo "Fatal flaw detected."
sleep 3  # Waits for 3 seconds
echo "Attempting possible workaround."
sleep 3  # Waits for 3 seconds
echo "Consulting alien life forms."
sleep 3  # Waits for 3 seconds
echo "Solution identified."
sleep 3  # Waits for 3 seconds
echo "Installation complete. Environment setup successfully completed."