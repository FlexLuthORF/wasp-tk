#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import shutil  # For copying files

def process_files(input_csv, original_output_csv, modified_output_csv):
    # Copy the original output CSV file to the modified output path
    shutil.copyfile(original_output_csv, modified_output_csv)
    
    # Load the copied file (now the one to be modified) and rename "genotyper_gene" to "gene"
    output_df = pd.read_csv(modified_output_csv)
    if 'genotyper_gene' in output_df.columns:
        output_df.rename(columns={'genotyper_gene': 'gene'}, inplace=True)

    # Load the input CSV
    input_df = pd.read_csv(input_csv)

    # Iterate over each row in the input DataFrame
    for _, input_row in input_df.iterrows():
        # Check for NA in 'allele_sequence_start' or 'allele_sequence_end' before proceeding
        if pd.notna(input_row['allele_sequence_start']) and pd.notna(input_row['allele_sequence_end']):
            # Match based on 'contig' and 'gene'
            matched_rows = output_df[(output_df['contig'] == input_row['contig']) & (output_df['gene'] == input_row['gene'])]
            
            # Append 'allele_sequence_start' and 'allele_sequence_end' to the matched rows in the output DataFrame
            for idx in matched_rows.index:
                output_df.at[idx, 'allele_sequence_start'] = input_row['allele_sequence_start']
                output_df.at[idx, 'allele_sequence_end'] = input_row['allele_sequence_end']

    # Save the modified DataFrame back to the new output CSV
    output_df.to_csv(modified_output_csv, index=False)

# Command-line arguments for input file path, original output file path, and modified output file path
input_csv_path = sys.argv[1]
original_output_csv_path = sys.argv[2]
modified_output_csv_path = sys.argv[3]

# Call the function with command-line arguments
process_files(input_csv_path, original_output_csv_path, modified_output_csv_path)
