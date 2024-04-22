#!/bin/bash

set -e -x

file=$(ls *fofn* | head -n 1)
user=$(whoami)
dir=./run_hifiasm
cat $file | while read sample data_path
do
    mkdir -p ${dir}/${sample}/merged_bam
    mkdir -p ${dir}/${sample}/merged_bam/alg_asm20_to_ref
    mkdir -p ${dir}/${sample}/merged_bam/alg_asm20_to_ref_with_secondarySeq
    outdir=${dir}/${sample}/merged_bam/alg_asm20_to_ref_with_secondarySeq

    sbatch --time=72:00:00 -p compute -o ./jobs/${sample}_process_job.txt --wrap="source /home/zmvanw01/git_repos/wasp/hifi-mapping/new_merge.sh; merge_and_rmdup ${sample}; align_and_process ${sample}"

    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 25 ]
    do
        sleep 1s
        count=`squeue | grep $user | wc -l`
    done
done