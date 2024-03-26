import argparse
import subprocess
import time

# Set up argument parsing
parser = argparse.ArgumentParser(description='Process sample IDs.')
parser.add_argument('sample_ids_file', type=str, help='Path to the file containing sample IDs')
args = parser.parse_args()

# Read sample IDs from the file
sample_ids = []
with open(args.sample_ids_file, 'r') as file:
    for line in file:
        # Assuming the file is tab-delimited; adjust as needed
        sample_id = line.split()[0]
        sample_ids.append(sample_id)

# Path to the script
script_path = '/home/zmvanw01/git_repos/swrm_scripts/zvw/repseq_personal-igblast.py'

# Maximum number of concurrent jobs
max_jobs = 12

# Function to get the current number of running jobs
def get_running_jobs_count():
    result = subprocess.run(['squeue', '-u', 'zmvanw01'], capture_output=True, text=True)
    return len(result.stdout.splitlines()) - 1  # Adjust as needed

# Iterate over sample IDs and submit jobs
for sample_id in sample_ids:
    # Wait if the maximum number of jobs is reached
    while get_running_jobs_count() >= max_jobs:
        time.sleep(60)  # Wait for 60 seconds before checking again

    # Construct the command to be wrapped
    command_to_wrap = f"python -u {script_path} {sample_id}"

    # Submit the job with --wrap
    log_file = f"./logs/{sample_id}.log"
    subprocess.run(["sbatch", 
                    "--job-name=repseq", 
                    "-o", log_file, 
                    "--wrap", command_to_wrap])

    print(f"Submitted job for sample ID: {sample_id}")
