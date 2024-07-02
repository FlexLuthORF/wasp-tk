#!/usr/bin/env python3

import pandas as pd
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
        # Define the start and end column names based on the 'gene' column value
        if any(substring in input_row['gene'] for substring in ['IGKV', 'IGLV', 'IGHV']):
            start_col, end_col = 'V-REGION_start', 'V-REGION_end'
        elif any(substring in input_row['gene'] for substring in ['IGKJ', 'IGLJ', 'IGHJ']):
            start_col, end_col = 'J-REGION_start', 'J-REGION_end'
        elif 'IGHD' in input_row['gene']:
            start_col, end_col = 'D-REGION_start', 'D-REGION_end'
        elif any(substring in input_row['gene'] for substring in ['IGKC', 'IGLC']):
            start_col, end_col = 'allele_sequence_start', 'allele_sequence_end'
        else:
            # Skip this row if none of the conditions are met
            continue

        # Check for NA in the determined start or end columns before proceeding
        if pd.notna(input_row[start_col]) and pd.notna(input_row[end_col]):
            # Match based on 'contig' and 'gene'
            matched_rows = output_df[(output_df['contig'] == input_row['contig']) & (output_df['gene'] == input_row['gene'])]
            
            # Append the determined start and end values to the matched rows in the output DataFrame
            for idx in matched_rows.index:
                output_df.at[idx, 'REGION_start'] = input_row[start_col]
                output_df.at[idx, 'REGION_end'] = input_row[end_col]

    # Replace all commas in the 'notes' column with semicolons
    if 'notes' in output_df.columns:
        output_df['notes'] = output_df['notes'].str.replace(',', ';', regex=False)

    # Save the modified DataFrame back to the new output CSV
    output_df.to_csv(modified_output_csv, index=False)

# Command-line arguments for input file path, original output file path, and modified output file path
input_csv_path = sys.argv[1]
original_output_csv_path = sys.argv[2]
modified_output_csv_path = sys.argv[3]

# Call the function with command-line arguments
process_files(input_csv_path, original_output_csv_path, modified_output_csv_path)
