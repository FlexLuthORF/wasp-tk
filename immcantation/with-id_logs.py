import os
import pandas as pd
import sys

# Check for the input argument
if len(sys.argv) < 2:
    print("Usage: script.py <input_dirs_file>")
    sys.exit(1)

input_dirs_file = sys.argv[1]

# Read the first column (directory names) from the file
input_dirs = []
with open(input_dirs_file, 'r') as file:
    for line in file:
        dir_name = line.split()[0]  # Split the line and take the first element
        input_dirs.append(dir_name)

output_dir = "./master/with-id_logs"
os.makedirs(output_dir, exist_ok=True)

for dir_name in input_dirs:
    log_dir = os.path.join(dir_name, "logs")
    if os.path.exists(log_dir):
        for filename in os.listdir(log_dir):
            if filename.endswith(".tab"):
                file_path = os.path.join(log_dir, filename)
                df = pd.read_csv(file_path, sep='\t', low_memory=False)
                df['SampleID'] = dir_name
                output_file = os.path.join(output_dir, filename)
                if os.path.exists(output_file):
                    df.to_csv(output_file, mode='a', header=False, sep='\t', index=False)
                else:
                    df.to_csv(output_file, mode='w', header=True, sep='\t', index=False)
