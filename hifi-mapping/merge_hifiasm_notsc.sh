#!/bin/bash

set -e -x

reffn=/home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta

dir=./run_hifiasm

# Read the sample names from the first column of the input file
while read sample; do
    file1="${dir}/${sample}/hifiasm/asm.bp.hap1.p_ctg_to_ref.sorted.bam"
    file2="${dir}/${sample}/hifiasm/asm.bp.hap2.p_ctg_to_ref.sorted.bam"

    if [ -f "$file1" ] && [ -f "$file2" ]; then
        mkdir -p ${dir}/${sample}/merged_bam_notsc
        mkdir -p ${dir}/${sample}/merged_bam/alg_asm20_notsc_to_ref
        outdir=${dir}/${sample}/merged_bam/alg_asm20_notsc_to_ref

        mkdir -p ${outdir}/jobs

        sbatch -o ${outdir}/jobs/align_job_notsc.txt --wrap="minimap2 -x asm20 -t 10 -L -a ${reffn} \
              ${dir}/${sample}/merged_bam_notsc/merged_all_reads.rmdup.fasta > ${outdir}/${sample}.sam
              wait
              samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
              samtools sort -@ 10 ${outdir}/${sample}.bam -o ${outdir}/${sample}.sorted.bam
              samtools index ${outdir}/${sample}.sorted.bam"

    else
        echo "One or both files for sample ${sample} do not exist."
    fi

    user=$(whoami)
    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 25 ]; do
        sleep 1s
        count=`squeue | grep $user | wc -l`
    done
done < <(cut -f1 "$1")