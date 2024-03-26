import subprocess
import pandas as pd
import os
import argparse
import glob

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('sampleID', type=str, help='Sample ID')
args = parser.parse_args()
prestodir = "/home/watsonlab/repSeq/CW50_NS_pool2/presto"
sampleID = args.sampleID
changeo_folder = f"{prestodir}/{sampleID}/changeo"
sample_seq_file = f"{prestodir}/{sampleID}/S5-final_total_collapse-unique_atleast-2_reheader.fasta"

for chain in ('IGK','IGL'):
    igblast_output = os.path.join(changeo_folder, f"igblast_output_{chain}.fmt7")
    if not os.path.exists(igblast_output):
        database_dir = f'{prestodir}/{sampleID}/alleles/personal-ref/{chain}/database'
        fasta_dir = f'{prestodir}/{sampleID}/alleles/personal-ref/{chain}/fasta'
        igdata_path = f'{prestodir}/{sampleID}/alleles/personal-ref/{chain}'
        os.makedirs(igdata_path, exist_ok=True)
        os.makedirs(database_dir, exist_ok=True)
        os.makedirs(fasta_dir, exist_ok=True)

        # Copying files
        subprocess.run(["cp", "-r", "/home/zmvanw01/share/igblast/internal_data", igdata_path])
        subprocess.run(["cp", "-r", "/home/zmvanw01/share/igblast/optional_file", igdata_path])
        d_allele_files = glob.glob('/home/zmvanw01/share/igblast/database/imgt_human_ig_d.*')
        for file in d_allele_files:
            subprocess.run(["cp", file, database_dir])

        # Run fasta-from-annotations.py script
        subprocess.run(["python", "/home/zmvanw01/git_repos/swrm_scripts/zvw/fasta-from-annotations.py",
                        f"/home/egenge01/projects/CW50/{chain}_alleles/{sampleID}/annotations.csv", 
                        fasta_dir])
        
        # Setting up file paths
        v_alleles_file = f"{fasta_dir}/human_ig_V_alleles.fasta"
        d_alleles_file = f"{fasta_dir}/human_ig_D_alleles.fasta"
        j_alleles_file = f"{fasta_dir}/human_ig_J_alleles.fasta"
        

        # Create custom BLAST databases
        subprocess.run(["makeblastdb", "-in", v_alleles_file, "-dbtype", "nucl", "-parse_seqids", "-out", f"{database_dir}/human_ig_V"])
        subprocess.run(["makeblastdb", "-in", j_alleles_file, "-dbtype", "nucl", "-parse_seqids", "-out", f"{database_dir}/human_ig_J"])
        
        # Setup for igblastn
        
        os.makedirs(changeo_folder, exist_ok=True)
        

        # Custom environment for igblastn
        igblast_env = os.environ.copy()
        print("path: ", igdata_path)
        igblast_env["IGDATA"] = igdata_path
        print(igblast_env["IGDATA"])
        igblast_command = [
            "igblastn",
            "-germline_db_V", f"{database_dir}/human_ig_V",
            "-germline_db_D", f"{database_dir}/imgt_human_ig_d",
            "-germline_db_J", f"{database_dir}/human_ig_J",
            "-auxiliary_data", f"{igdata_path}/optional_file/human_gl.aux",
            "-query", sample_seq_file,
            "-outfmt", "19",
            "-out", igblast_output,
            "-organism", "human",
            "-ig_seqtype", "Ig",
            "-domain_system", "imgt"
        ]
        subprocess.run(igblast_command, env=igblast_env)
    else:
        print("skipping ", chain)
# Handling IGH chain
igdata_path = '/home/zmvanw01/share/IGH'
print("path: ", igdata_path)
igblast_env = os.environ.copy()
igblast_env["IGDATA"] = igdata_path
print(igblast_env["IGDATA"])
igblast_output = os.path.join(changeo_folder, "igblast_output_H.fmt7")
if not os.path.exists(igblast_output):
    igblast_command = [
    "AssignGenes.py", "igblast",
    "-s", sample_seq_file,
    "-b", "/home/zmvanw01/share/igblast",
    "--format", "airr",
    "-o", igblast_output,
    "--organism", "human",
    "--loci", "ig"
    ]
    subprocess.run(igblast_command)