#!/bin/bash

# Check if at least one argument was provided (the path to fofn.tsv)
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 path/to/your/fofn.tsv [dir1 dir2 ...]"
    exit 1
fi

# Use the first argument as the path to your fofn.tsv file
fofn_path="$1"

# Shift the arguments so $@ does not include the first one (path to fofn.tsv)
shift

# Base directory where the CSV files are located
base_dir=$(pwd)

# Default array of directories to loop through
dirs=("igh" "igk" "igl" "ighc")

# Check if any additional arguments were provided
if [ "$#" -gt 0 ]; then
    # If yes, override the dirs array with the provided arguments
    dirs=("$@")
fi

# Loop through each directory
for dir in "${dirs[@]}"; do
    # Set the command options based on the directory
    case $dir in
        igh)
            options=("IGH" "igh" "+-")
            ;;
        igk)
            options=("IGK" "chr2" "+-")
            ;;
        igl)
            options=("IGL" "chr22" "+")
            ;;
        ighc)
            options=("IGH" "ighc" "+-")
            ;;                      
        *)
            echo "Warning: Unsupported directory $dir. Skipping..."
            continue
            ;;
    esac

    # Directory for current iteration
    csv_dir="${base_dir}/${dir}"

    # Loop through the fofn.tsv file
    while IFS=$'\t' read -r sampleID _; do
        # Use eval to correctly handle the inclusion of quotes in options
        echo python /home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/import_from_assemblies.py "${options[0]}" "${options[1]}" "${options[2]}" "${csv_dir}/${sampleID}_genes.csv" /home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta /home/zmvanw01/git_repos/immune_receptor_genomics/current /home/zmvanw01/git_repos/VDJbase-Genomics/ref/ "${csv_dir}/${sampleID}_genes_imported.csv"
        eval python /home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/import_from_assemblies.py "${options[0]}" "${options[1]}" "${options[2]}" "${csv_dir}/${sampleID}_genes.csv" /home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta /home/zmvanw01/git_repos/immune_receptor_genomics/current /home/zmvanw01/git_repos/VDJbase-Genomics/ref/ "${csv_dir}/${sampleID}_genes_imported.csv"
    done < "$fofn_path"
done
