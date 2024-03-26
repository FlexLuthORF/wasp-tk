import pandas as pd
import os
import sys

def merge_files(sample_file):
    with open(sample_file, 'r') as file:
        for line in file:
            sampleID = line.split()[0].strip()
            changeo_path = f'./presto/{sampleID}/changeo'
            merged_file_path = f'{changeo_path}/{sampleID}_merged-changeo.tsv'

            # Skip if merged file exists
            if os.path.exists(merged_file_path):
                continue

            file_paths = [
                f'{changeo_path}/igblast_output_H.fmt7',
                f'{changeo_path}/igblast_output_IGK.fmt7',
                f'{changeo_path}/igblast_output_IGL.fmt7'
            ]

            merged_data = pd.DataFrame()
            for file_path, locus in zip(file_paths, ['IGH', 'IGK', 'IGL']):
                if os.path.exists(file_path):
                    data = pd.read_csv(file_path, sep='\t')
                    merged_data = pd.concat([merged_data, data[data['locus'] == locus]])

            if not merged_data.empty:
                os.makedirs(changeo_path, exist_ok=True)
                merged_data.to_csv(merged_file_path, sep='\t', index=False)

if len(sys.argv) > 1:
    merge_files(sys.argv[1])
else:
    print("Please provide a sample file name.")
