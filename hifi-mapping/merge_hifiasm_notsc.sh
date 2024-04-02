#!/bin/bash
set -e -x

#reffn=/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/reference.fasta
reffn=../../immune_receptor_genomics/current/reference.fasta
dir=./run_hifiasm
#cat $PWD/138_samples.txt | while read sample

do
    #    if [ ! -s ${dir}/${sample}/merged_bam/alg_asm20_notsc_to_ref/${sample}.sorted.bam ]
    #    then
    file1="${dir}/${sample}/hifiasm/asm.bp.hap1.p_ctg_to_ref.sorted.bam"
    file2="${dir}/${sample}/hifiasm/asm.bp.hap2.p_ctg_to_ref.sorted.bam"
    
    if [ -f "$file1" ] && [ -f "$file2" ]; then
	mkdir -p ${dir}/${sample}/merged_bam_notsc
	mkdir -p ${dir}/${sample}/merged_bam/alg_asm20_notsc_to_ref
	outdir=${dir}/${sample}/merged_bam/alg_asm20_notsc_to_ref
	mkdir -p ${outdir}/jobs
	sbatch --time=72:00:00 -c 1 -p compute -o ${outdir}/jobs/align_job_notsc.txt --wrap="minimap2 -x asm20 -t 10 -L -a ${reffn} \
    	  ${dir}/${sample}/merged_bam_notsc/merged_all_reads.rmdup.fasta > ${outdir}/${sample}.sam
    	  wait
    	  samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam 
    	  samtools sort -@ 10 ${outdir}/${sample}.bam -o ${outdir}/${sample}.sorted.bam
    	  samtools index ${outdir}/${sample}.sorted.bam"
#    	 rm -f ${outdir}/${sample}.sam
#    	 rm -f ${outdir}/${sample}.bam
    else
        echo "One or both files for sample ${sample} do not exist."
    fi
    #    fi
    
    # 	sbatch --time=72:00:00 -c 1 -p compute -o ${outdir}/jobs/align_job_notsc.txt --wrap="
    #       source activate IGv2
    #       samtools merge -f ${dir}/${sample}/merged_bam_notsc/merged.bam ${dir}/${sample}/hifiasm/asm.bp.hap1.p_ctg_to_ref.sorted.bam \
    #       ${dir}/${sample}/hifiasm/asm.bp.hap2.p_ctg_to_ref.sorted.bam
    #        wait
    #        samtools sort ${dir}/${sample}/merged_bam_notsc/merged.bam -o ${dir}/${sample}/merged_bam_notsc/merged.sorted.bam
    #        wait
    #        samtools index ${dir}/${sample}/merged_bam_notsc/merged.sorted.bam
    #        wait
    #        samtools fasta --reference ${reffn} ${dir}/${sample}/merged_bam_notsc/merged.bam > ${dir}/${sample}/merged_bam_notsc/merged_all_reads.fasta
    #        wait
    #        conda deactivate
    #        source activate seqkit-env
    #        seqkit rmdup --by-seq ${dir}/${sample}/merged_bam_notsc/merged_all_reads.fasta \
    #        -o ${dir}/${sample}/merged_bam_notsc/merged_all_reads.rmdup.fasta"
    user=$(whoami)
    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 15 ]
    do
	sleep 1s
	count=`squeue | grep $user | wc -l`
    done
done
