import os
import argparse
import pandas as pd

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert CSV allele files to FASTA format.')
    parser.add_argument('--directory', 
                        type=str, 
                        default='.', 
                        help='Directory containing input CSV files (default: current directory)')
    
    # Define input and output files
    input_files = [
        "master_TRA_alleles.csv",
        "master_TRB_alleles.csv",
        "master_TRD_alleles.csv",
        "master_TRG_alleles.csv",
        "master_IGH_alleles.csv",
        "master_IGHC_alleles.csv",
        "master_IGK_alleles.csv",
        "master_IGL_alleles.csv"
    ]
    
    fasta_outputs = [
        "TRA_alleles.fasta",
        "TRB_alleles.fasta",
        "TRD_alleles.fasta",
        "TRG_alleles.fasta",
        "IGH_alleles.fasta",
        "IGHC_alleles.fasta",
        "IGK_alleles.fasta",
        "IGL_alleles.fasta"
    ]
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process each file
    for input_file, fasta_output in zip(input_files, fasta_outputs):
        
        input_path = os.path.join(args.directory, input_file)
        fasta_output_path = os.path.join(args.directory, fasta_output)
        
        # Skip if input file doesn't exist
        if not os.path.exists(input_path):
            print(f"Warning: {input_file} not found. Skipping.")
            continue
        
        # Read CSV file
        df = pd.read_csv(input_path)
        
        # Apply filtering
        is_c_gene = df["vdjbase_allele"].str[3] == "C"
        c_filter = df["Average_Coverage"] >= 50
        non_c_filter = (
            (df["Average_Coverage"] >= 30) &
            (df["Fully_Spanning_Reads_100%_Match"] >= 10)
        )
        combined_filter = (is_c_gene & c_filter) | (~is_c_gene & non_c_filter)
        filtered = df[combined_filter]
        
        # Remove duplicates
        unique_entries = filtered.drop_duplicates(subset=["vdjbase_allele", "gene_sequence"])
        
        # Write FASTA file
        with open(fasta_output_path, "w") as fasta:
            for _, row in unique_entries.iterrows():
                fasta.write(f">{row['vdjbase_allele']}\n{row['gene_sequence']}\n")
        
        print(f"Processed {input_file} -> {fasta_output}")

if __name__ == "__main__":
    main()