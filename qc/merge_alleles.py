import os
import pandas as pd

def create_master_csv_files(start_dir):
    alleles_dir = 'alleles'
    loci_dirs = ["IGH", "IGHC", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]
    missing_samples = {}
    # Create a dictionary to store DataFrames for each loci
    master_data = {locus: [] for locus in loci_dirs}

    # Traverse all sample directories under the start directory
    for sample_id in os.listdir(start_dir):
        sample_dir = os.path.join(start_dir, sample_id)

        if os.path.isdir(sample_dir):
            alleles_path = os.path.join(sample_dir, alleles_dir)

            if os.path.exists(alleles_path):
                for loci in loci_dirs:
                    files_found = False

                    # Check each file in the alleles directory
                    for file in os.listdir(alleles_path):
                        if loci in file and file.endswith('.csv'):
                            files_found = True
                            file_path = os.path.join(alleles_path, file)
                            df = pd.read_csv(file_path)
                            master_data[loci].append(df)

                    if not files_found:
                        missing_samples.setdefault(loci, []).append(sample_id)

    # After processing all samples, merge and write master CSVs for each loci
    for loci, dfs in master_data.items():
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # Drop duplicate columns if necessary
            combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]
            master_file_path = os.path.join(start_dir, f"master_{loci}_alleles.csv")
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
