#!/bin/bash

# Step 1: Download and install Miniconda

# Download the latest Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Make the installer script executable
chmod +x Miniconda3-latest-Linux-x86_64.sh

# Run the Miniconda installation script
./Miniconda3-latest-Linux-x86_64.sh

# Step 2: Install dependencies using conda

# Install dependencies from conda-forge and bioconda channels
conda install -c conda-forge -c bioconda \
    fastqc \
    multiqc \
    fastp \
    bwa \
    samtools \
    bcftools

# Step 3: Exit the script with a success message
echo "Installation complete. All dependencies have been installed."

