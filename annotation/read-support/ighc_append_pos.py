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

    # Define the columns to be extracted and appended for start and end
    exon_columns = [
        'C-EXON_1_start', 'C-EXON_1_end',
        'C-EXON_2_start', 'C-EXON_2_end',
        'C-EXON_3_start', 'C-EXON_3_end',
        'C-EXON_4_start', 'C-EXON_4_end',
        'C-EXON_5_start', 'C-EXON_5_end',
        'C-EXON_6_start', 'C-EXON_6_end',
        'C-EXON_7_start', 'C-EXON_7_end',
        'C-EXON_8_start', 'C-EXON_8_end',
        'C-EXON_9_start', 'C-EXON_9_end'
    ]

    # Process each row in the input DataFrame
    for _, input_row in input_df.iterrows():
        # Match based on 'contig' and 'gene'
        matched_rows = output_df[(output_df['contig'] == input_row['contig']) & (output_df['gene'] == input_row['gene'])]

        # Initialize variables to store the minimum start and maximum end
        min_start = float('inf')
        max_end = 0

        # Collect all start and end values, find the minimum and maximum
        for col in exon_columns:
            if col.endswith("_start") and col in input_row and pd.notna(input_row[col]):
                start_value = input_row[col]
                if start_value > 0:  # Ensure the value is non-zero
                    min_start = min(min_start, start_value)

            if col.endswith("_end") and col in input_row and pd.notna(input_row[col]):
                max_end = max(max_end, input_row[col])

            # Append exon columns to the matched rows in the output DataFrame
            if col in input_row and pd.notna(input_row[col]):
                for idx in matched_rows.index:
                    output_df.at[idx, col] = input_row[col]

        # Assign the computed values to matched rows
        for idx in matched_rows.index:
            if min_start < float('inf'):
                output_df.at[idx, 'allele_sequence_start'] = min_start
            if max_end > 0:
                output_df.at[idx, 'allele_sequence_end'] = max_end

    # Save the modified DataFrame back to the new output CSV
    output_df.to_csv(modified_output_csv, index=False)

# Command-line arguments for input file path, original output CSV file path, and modified output CSV file path
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py input_csv original_output_csv modified_output_csv")
        sys.exit(1)

    input_csv_path = sys.argv[1]
    original_output_csv_path = sys.argv[2]
    modified_output_csv_path = sys.argv[3]

    process_files(input_csv_path, original_output_csv_path, modified_output_csv_path)