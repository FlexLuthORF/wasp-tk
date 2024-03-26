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

    # Define the columns to be extracted and appended
    exon_columns = [
        'C-EXON_4_start', 'C-EXON_4_end',
        'C-EXON_3_start', 'C-EXON_3_end',
        'C-EXON_2_start', 'C-EXON_2_end',
        'C-EXON_1_start', 'C-EXON_1_end',
        'C-EXON_6_start', 'C-EXON_6_end',
        'C-EXON_5_start', 'C-EXON_5_end',
        'C-EXON_9_start', 'C-EXON_9_end',
        'C-EXON_8_start', 'C-EXON_8_end',
        'C-EXON_7_start', 'C-EXON_7_end'
    ]

    # Iterate over each row in the input DataFrame
    for _, input_row in input_df.iterrows():
        # Match based on 'contig' and 'gene'
        matched_rows = output_df[(output_df['contig'] == input_row['contig']) & (output_df['gene'] == input_row['gene'])]
        
        # Append exon columns to the matched rows in the output DataFrame
        for idx in matched_rows.index:
            for col in exon_columns:
                # Check if the column exists in the input row and is not NA before appending
                if col in input_row and pd.notna(input_row[col]):
                    output_df.at[idx, col] = input_row[col]

    # Save the modified DataFrame back to the new output CSV
    output_df.to_csv(modified_output_csv, index=False)

# Command-line arguments for input file path, original output CSV file path, and modified output CSV file path
input_csv_path = sys.argv[1]
original_output_csv_path = sys.argv[2]
modified_output_csv_path = sys.argv[3]

# Call the function with command-line arguments
process_files(input_csv_path, original_output_csv_path, modified_output_csv_path)
