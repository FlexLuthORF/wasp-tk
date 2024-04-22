#!/bin/bash
set -e -x

#reffn=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/reference.fasta
reffn=/home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta
dir=./run_hifiasm
#cat /home/egenge01/projects/CW50/filtered_2023-10-25_fofn.txt | while read sample data_path 
#do
#cat batch2_samples_fofn.txt | grep '220880402C' | while read sample data_path
#cat Seq_CW2_7_CW40_23_CW50_27-32_paths.fofn.txt | while read sample data_path
file=$(ls *fofn* | head -n 1)
user=$(whoami)
cat $file | while read sample data_path
do
#    if [ ! -s ${dir}/${sample}/merged_bam/alg_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ]
#    then
    mkdir -p ${dir}/${sample}/merged_bam
    mkdir -p ${dir}/${sample}/merged_bam/alg_asm20_to_ref
    mkdir -p ${dir}/${sample}/merged_bam/alg_asm20_to_ref_with_secondarySeq
    outdir=${dir}/${sample}/merged_bam/alg_asm20_to_ref_with_secondarySeq
    mkdir -p ${outdir}/jobs
    sbatch --time=72:00:00 -p compute -o ${outdir}/jobs/align_job2.txt --wrap="minimap2 -x asm20 --secondary-seq -t 10 -L -a ${reffn} ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta > ${outdir}/${sample}.sam;
			 wait
			 samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam; 
			 samtools sort -@ 10 ${outdir}/${sample}.bam -o ${outdir}/${sample}.sorted.bam;
			 samtools index ${outdir}/${sample}.sorted.bam
			 #rm -f ${outdir}/${sample}.sam;
			 #rm -f ${outdir}/${sample}.bam
"
	#	fi
    #conda activate seqkit-env
    sbatch --time=72:00:00 -p compute -o ${dir}/${sample}/merged_bam/rmdup_job.txt --wrap="samtools merge -f ${dir}/${sample}/merged_bam/merged.bam ${dir}/${sample}/break_at_soft_clip/2/1_asm20_hifi_asm_to_ref.sorted.bam ${dir}/${sample}/break_at_soft_clip/2/2_asm20_hifi_asm_to_ref.sorted.bam;
wait
samtools sort ${dir}/${sample}/merged_bam/merged.bam -o ${dir}/${sample}/merged_bam/merged.sorted.bam;
wait
samtools index ${dir}/${sample}/merged_bam/merged.sorted.bam;
samtools fasta --reference ${reffn} ${dir}/${sample}/merged_bam/merged.bam > ${dir}/${sample}/merged_bam/merged_all_reads.fasta
seqkit rmdup --by-seq ${dir}/${sample}/merged_bam/merged_all_reads.fasta -o ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta"

   samtools merge -f ${dir}/${sample}/merged_bam/merged.bam ${dir}/${sample}/break_at_soft_clip/2/1_asm20_hifi_asm_to_ref.sorted.bam ${dir}/${sample}/break_at_soft_clip/2/2_asm20_hifi_asm_to_ref.sorted.bam;
   wait
   samtools sort -@ 8 ${dir}/${sample}/merged_bam/merged.bam -o ${dir}/${sample}/merged_bam/merged.sorted.bam;
   wait
   samtools index ${dir}/${sample}/merged_bam/merged.sorted.bam;
   samtools fasta --reference ${reffn} ${dir}/${sample}/merged_bam/merged.sorted.bam > ${dir}/${sample}/merged_bam/merged_all_reads.fasta
   seqkit rmdup --by-seq ${dir}/${sample}/merged_bam/merged_all_reads.fasta -o ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta
       fi
   count=`squeue | grep egenge01 | wc -l`
   while [ ${count} -gt 40 ]
   do
	sleep 1s
	count=`squeue | grep egenge01 | wc -l`
   done
done
