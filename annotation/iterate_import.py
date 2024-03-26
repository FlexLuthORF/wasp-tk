import os
import subprocess
import sys
import time

# Check if an argument was provided
if len(sys.argv) != 2:
    print("Usage: python3 this_script.py path/to/your/fofn.tsv")
    sys.exit(1)

# Use the first argument as the path to your fofn.tsv file
fofn_path = sys.argv[1]

# Create necessary directories
sbatch_dir = "./slurm/sbatch"
logs_dir = "./slurm/logs"
os.makedirs(sbatch_dir, exist_ok=True)
os.makedirs(logs_dir, exist_ok=True)

# Dictionary of directories and their corresponding options
dirs = {
    #"igh": ("IGH", "igh", "+-"),
    #"igk": ("IGK", "chr2", "+-"),
    #"igl": ("IGL", "chr22", "+"),
    "ighc": ("IGH", "ighc", "+-")
}

def check_jobs(user):
    """Check the number of running jobs for a specific user."""
    result = subprocess.run(['squeue', '-u', user, '-h'], capture_output=True, text=True)
    return len(result.stdout.splitlines())

def create_sbatch_command(sampleID):
    """Constructs the SBATCH command for a given sampleID."""
    commands = []
    for dir, options in dirs.items():
        commands.append(f"python /home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/import_from_assemblies.py {options[0]} {options[1]} {options[2]} {dir}/{sampleID}_genes.csv /home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta /home/zmvanw01/git_repos/immune_receptor_genomics/current /home/zmvanw01/git_repos/VDJbase-Genomics/ref/ {dir}/{sampleID}_genes_imported.csv")

    sbatch_script_content = f"""#!/bin/bash
#SBATCH --job-name={sampleID}_job
#SBATCH --output={logs_dir}/{sampleID}.out
#SBATCH --error={logs_dir}/{sampleID}.err
#SBATCH --nodes=1

{' && '.join(commands)}
"""
    sbatch_script_path = f"{sbatch_dir}/{sampleID}.sbatch"
    with open(sbatch_script_path, 'w') as f:
        f.write(sbatch_script_content)
    
    return sbatch_script_path

# Main loop
with open(fofn_path, 'r') as f:
    for line in f:
        sampleID = line.split('\t')[0].strip()

        # Ensure we don't exceed the job limit
        while check_jobs('zmvanw01') >= 12:
            print("Job limit reached. Waiting...")
            time.sleep(10)  # Wait for 1 minute before checking again

        # Create and submit the SBATCH script
        sbatch_script_path = create_sbatch_command(sampleID)
        subprocess.run(['sbatch', sbatch_script_path])
        print(f"Submitted: {sbatch_script_path}")
