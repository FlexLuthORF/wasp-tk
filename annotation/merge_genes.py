import pandas as pd
import os

# Directory names
directories = ['igh', 'igk', 'igl']
output_dir = 'merge'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to get unique sample IDs from one directory (assuming consistent naming across directories)
def get_sample_ids_from_dir(directory):
    return [filename.split('_')[0] for filename in os.listdir(directory) if filename.endswith('_genes.csv')]

# Assuming 'igh' directory has all sample IDs
sample_ids = get_sample_ids_from_dir('igh')

# Process and merge CSV files for each sample ID
for sample_id in sample_ids:
    dfs = []  # List to hold DataFrames for merging
    
    # Read and append DataFrame from each directory
    for directory in directories:
        filepath = os.path.join(directory, f"{sample_id}_genes.csv")
        if os.path.exists(filepath):
            df = pd.read_csv(filepath)
            dfs.append(df)
    
    # Merge DataFrames, retaining all columns
    merged_df = dfs[0]  # Start with the first DataFrame (assumed to be 'igh' with all columns)
    for df in dfs[1:]:  # Merge with remaining DataFrames
        merged_df = pd.merge(merged_df, df, how='outer', on='gene')
    
    # Save the merged DataFrame
    output_filepath = os.path.join(output_dir, f"{sample_id}_all-genes.csv")
    merged_df.to_csv(output_filepath, index=False)

print("Merging complete. Files saved in the 'merge' directory.")
