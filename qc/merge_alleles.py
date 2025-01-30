import os
import pandas as pd

def create_master_csv_files(start_dir):
    alleles_dir = 'alleles'
    loci_dirs = ["IGH", "IGHC", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]
    missing_samples = {}

    # Traverse all sampleId directories under the start directory
    for sample_id in os.listdir(start_dir):
        sample_dir = os.path.join(start_dir, sample_id)
        
        if os.path.isdir(sample_dir):  # Only proceed if it's a directory
            alleles_path = os.path.join(sample_dir, alleles_dir)

            # If the alleles directory exists
            if os.path.exists(alleles_path):
                for loci in loci_dirs:
                    master_file_path = os.path.join(start_dir, f"master_{loci}_alleles.csv")
                    combined_df = None
                    files_found = False

                    # Traverse all files in the alleles directory of this sample
                    for file in os.listdir(alleles_path):
                        if loci in file and file.endswith('.csv'):
                            files_found = True
                            file_path = os.path.join(alleles_path, file)
                            df = pd.read_csv(file_path)

                            # Combine the CSVs
                            if combined_df is None:
                                combined_df = df
                            else:
                                combined_df = pd.concat([combined_df, df], ignore_index=True)

                    # If no files were found for this loci in this sample, track it
                    if not files_found:
                        if loci not in missing_samples:
                            missing_samples[loci] = []
                        missing_samples[loci].append(sample_id)

                    # If we have combined data for the current loci, save it
                    if combined_df is not None:
                        # Drop duplicate columns
                        combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]
                        combined_df.to_csv(master_file_path, index=False)

    # Output missing sample IDs to a missing.txt file
    with open(os.path.join(start_dir, "missing.txt"), "w") as f:
        for loci, samples in missing_samples.items():
            f.write(f"Missing files for {loci} in the following samples:\n")
            for sample in samples:
                f.write(f"{sample}\n")
            f.write("\n")

if __name__ == "__main__":
    start_dir = os.getcwd()
    create_master_csv_files(start_dir)
