#!/bin/bash

set -e -x

user=$(whoami)

input_file=$(ls *fofn* | head -n 1)

while read -r sample; do
    bam_file="/$PWD/run_hifiasm/$sample/merged_bam/alg_asm20_to_ref_with_secondarySeq/$sample.sorted.bam"
    
    if [ -f "$bam_file" ]; then
        outdir="$PWD/geno_analysis/${sample}"
        mkdir -p "$outdir"
        
        sbatch --time=8:00:00 -p compute -o "${outdir}/job.txt" --wrap="sh /home/zmvanw01/git_repos/wasp/annotation/get_vcf/final_vcf.sh $sample"
        
        count=$(squeue | grep $user | wc -l)
        
        while [ "$count" -gt 9 ]; do
            sleep 1s
            count=$(squeue | grep $user | wc -l)
        done
    else
        echo "BAM file not found for sample: $sample"
    fi
done < <(cut -f1 "$input_file")