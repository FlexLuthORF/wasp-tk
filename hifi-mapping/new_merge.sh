#!/bin/bash

set -e -x

reffn=/home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta
#root = $PWD
dir=$PWD/run_hifiasm

function align_and_process {
    sample=$1
    outdir=${dir}/${sample}/merged_bam/alg_asm20_to_ref_with_secondarySeq
    mkdir -p $outdir
    minimap2 -x asm20 --secondary-seq -t 12 -L -a ${reffn} ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta > ${outdir}/${sample}.sam
    samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
    samtools sort -@ 10 ${outdir}/${sample}.bam -o ${outdir}/${sample}.sorted.bam
    samtools index ${outdir}/${sample}.sorted.bam
}

function merge_and_rmdup {
    sample=$1
    mkdir -p ${dir}/${sample}/merged_bam

    samtools merge -f ${dir}/${sample}/merged_bam/merged.bam ${dir}/${sample}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam ${dir}/${sample}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam
    samtools sort -@ 12 ${dir}/${sample}/merged_bam/merged.bam -o ${dir}/${sample}/merged_bam/merged.sorted.bam
    samtools index ${dir}/${sample}/merged_bam/merged.sorted.bam
    samtools fasta --reference ${reffn} ${dir}/${sample}/merged_bam/merged.sorted.bam > ${dir}/${sample}/merged_bam/merged_all_reads.fasta
    seqkit rmdup --by-seq ${dir}/${sample}/merged_bam/merged_all_reads.fasta -o ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta
}