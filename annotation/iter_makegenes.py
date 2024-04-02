import os
import argparse
import subprocess

def create_slurm_script(job_id, samples):
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name=process_samples_{job_id}
#SBATCH --output=process_samples_{job_id}_%j.out
#SBATCH --error=process_samples_{job_id}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12


"""
    for sample_id, bam_path in samples:
        slurm_script += f"python /home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/process_samples.py {sample_id} {bam_path}\n"

    return slurm_script

def main(fofn_file):
    with open(fofn_file, 'r') as f:
        samples = [line.strip().split('\t') for line in f]

    num_samples = len(samples)
    samples_per_job = 3
    num_jobs = (num_samples + samples_per_job - 1) // samples_per_job

    for job_id in range(num_jobs):
        start_index = job_id * samples_per_job
        end_index = min(start_index + samples_per_job, num_samples)
        job_samples = samples[start_index:end_index]

        slurm_script = create_slurm_script(job_id, job_samples)

        with open(f"process_samples_{job_id}.slurm", 'w') as f:
            f.write(slurm_script)

        subprocess.run(f"sbatch process_samples_{job_id}.slurm", shell=True, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Submit Slurm jobs to process samples using make_gene_file.py and import_from_assemblies.py')
    parser.add_argument('fofn_file', help='path to the file of filenames (fofn) containing sample IDs and BAM paths')
    args = parser.parse_args()

    main(args.fofn_file)