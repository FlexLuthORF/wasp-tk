#!/bin/bash

# Function to download reference files
fetch_reference_files() {
    echo "Downloading reference.fasta..."
    wget -O ./data/reference.fasta http://immunogenomics.louisville.edu/immune_receptor_genomics/current/reference.fasta

    echo "Cloning repository for BED files..."
    git clone git@github.com:Watson-IG/immune_receptor_genomics.git

    echo "Copying the most recent data from 'current'..."
    cp -r immune_receptor_genomics/current ./data

    echo "Cleaning up cloned repository..."
    rm -rf immune_receptor_genomics
}

# Create directories if they don't exist
mkdir -p ./container
mkdir -p ./data

# Fetch the Singularity container file
echo "Downloading hifi_container.sif..."
wget -O ./container/hifi_container.sif http://immunogenomics.louisville.edu/wasp/hifi_container.sif

# Check if --ref is passed as an argument
if [[ "$1" == "--ref" ]]; then
    fetch_reference_files
fi

echo "Setup completed."
