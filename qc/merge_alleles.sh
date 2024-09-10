#!/bin/bash

# Set the current working directory
PWD=$(pwd)

# Define the loci options
LOCI_OPTIONS=("IGH" "IGHC" "IGK" "IGL" "TRA" "TRB" "TRD" "TRG")

# Build a list of sampleIds based on subdirectories in PWD
sampleIds=($(ls -d */ | sed 's#/##'))

# Loop over each loci option
for LOCI in "${LOCI_OPTIONS[@]}"; do
    # Prepare the master file for the current loci
    master_file="${PWD}/master_${LOCI}.csv"
    touch "$master_file"  # Create the master file if it doesn't exist
    first_file=1  # Flag to check if it's the first file

    # Process each sampleId
    for sampleId in "${sampleIds[@]}"; do
        file_path="${PWD}/${sampleId}/*/read_support/${sampleId}/imported_genes/${LOCI}/${sampleId}_make_gene_file_imported_with_read_support.csv"

        # Check if the file exists and append to the master file; otherwise, log to missing.csv
        if [ -f "$file_path" ]; then
            if [ $first_file -eq 1 ]; then
                cat "$file_path" > "$master_file"
                first_file=0  # Set the first file flag to 0 after the first file is processed
            else
                tail -n +2 "$file_path" >> "$master_file"  # Append without the header
            fi
        else
            echo "$sampleId missing $file_path" >> "${PWD}/missing.csv"
        fi
    done
done