import os
import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Process .tab files from sample ID directories.')
parser.add_argument('sample_id_file', type=str, help='Path to the file containing sample IDs.')
args = parser.parse_args()

# Read sample IDs from the provided file (tab-separated .txt file)
sample_ids_df = pd.read_csv(args.sample_id_file, sep='\t', header=None)
input_dirs = sample_ids_df.iloc[:, 0].tolist()

# Define the output directory
output_dir = "./master/with-id_logs"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Process each .tab file
for dir_name in input_dirs:
    log_dir = os.path.join("./presto", dir_name, "logs")
    if os.path.exists(log_dir):
        for filename in os.listdir(log_dir):
            if filename.endswith(".tab"):
                # Read the .tab file
                file_path = os.path.join(log_dir, filename)
                df = pd.read_csv(file_path, sep='\t', low_memory=False)

                # Add the sample ID as a new column
                df['SampleID'] = dir_name

                # Define the output file path
                output_file = os.path.join(output_dir, filename)

                # Concatenate or create the file with headers only for the first file
                if os.path.exists(output_file):
                    df.to_csv(output_file, mode='a', header=False, sep='\t', index=False)
                else:
                    df.to_csv(output_file, mode='w', header=True, sep='\t', index=False)
