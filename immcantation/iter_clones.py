import subprocess
import time
import sys
import os

# Check if the fofn file is provided as a command-line argument
if len(sys.argv) < 2:
    print("Please provide the fofn file as a command-line argument.")
    sys.exit(1)

fofn_file = sys.argv[1]

# Read sample names from the fofn file
with open(fofn_file, "r") as fofn:
    samples = [line.strip() for line in fofn]

# Get the current user using whoami
current_user = subprocess.check_output("whoami", shell=True, text=True).strip()

# Iterate over each sample
for sample in samples:
    # Check the number of running jobs for the current user
    while True:
        squeue_output = subprocess.check_output(f"squeue -u {current_user}", shell=True, text=True)
        num_running_jobs = len(squeue_output.strip().split("\n")) - 1

        if num_running_jobs < 15:
            break
        else:
            print(f"Number of running jobs exceeds 15. Waiting for 60 seconds...")
            time.sleep(60)  # Wait for 60 seconds before checking again

    # Define the loci to process
    loci = ["H", "IGK", "IGL"]

    # Rename the files for each locus
    for locus in loci:
        file_path = f"./presto/changeo/{sample}/igblast_output_{locus}.fmt7"

        # Check if the file already ends with ".tsv"
        if not file_path.endswith(".tsv"):
            new_file_path = f"{file_path}.tsv"
            if os.path.exists(file_path):
                subprocess.run(f"mv {file_path} {new_file_path}", shell=True)

    # Construct the Rscript command
    rscript_command = f"Rscript /home/zmvanw01/git_repos/wasp/immcantation/new_clones.R {sample}"

    # Submit a Slurm job using sbatch
    slurm_command = f"sbatch --job-name={sample}_clones --output={sample}_clones.out --error={sample}_clones.err --wrap='{rscript_command}'"
    subprocess.run(slurm_command, shell=True)