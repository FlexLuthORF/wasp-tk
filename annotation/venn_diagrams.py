import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import glob
import os

# Variables to modify
hifiasm_path = "/home/zmvanw01/projects/12-sample-comparison/hifiasm-path/annotations/igh"
ig_path = "/home/zmvanw01/projects/12-sample-comparison/ig-path/annotations/igh"
plots_dir = "plots"  # Directory to save plots

# Ensure plots directory exists
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Build sample list by globbing
file_pattern = os.path.join(hifiasm_path, "*_genes_imported.csv")
file_paths = glob.glob(file_pattern)
sample_ids = list(set([os.path.basename(fp).split('_genes_imported.csv')[0] for fp in file_paths]))

def compare_vdjbase_allele(hifiasm_file, ig_file):
    # Load CSV files
    df_hifiasm = pd.read_csv(hifiasm_file)
    df_ig = pd.read_csv(ig_file)
    
    # Extract 'vdjbase_allele' values
    alleles_hifiasm = set(df_hifiasm['vdjbase_allele'].dropna())
    alleles_ig = set(df_ig['vdjbase_allele'].dropna())
    
    return alleles_hifiasm, alleles_ig

def plot_venn(alleles_hifiasm, alleles_ig, sample_id):
    plt.figure(figsize=(8, 8))
    venn2([alleles_hifiasm, alleles_ig], set_labels=('hifiasm', 'ig'))
    plt.title(f"Venn Diagram for {sample_id}")
    # Save the plot to the specified plots directory
    plt.savefig(f"{plots_dir}/{sample_id}_venn_diagram.png")
    plt.close()  # Close the plot to free memory

for sample_id in sample_ids:
    hifiasm_file = f"{hifiasm_path}/{sample_id}_genes_imported.csv"
    ig_file = f"{ig_path}/{sample_id}_genes_imported.csv"
    
    alleles_hifiasm, alleles_ig = compare_vdjbase_allele(hifiasm_file, ig_file)
    plot_venn(alleles_hifiasm, alleles_ig, sample_id)
