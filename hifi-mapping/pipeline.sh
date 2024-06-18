#!/bin/bash
set -e -x

outdir=$1
ccs=$2
threads=$3
sample=$4

function align_with_minimap2 {
    fasta=$1
    prefix=$2
    reffn=$3
    threads=$4
    minimap2 \
	-t ${threads} --secondary=yes -L -a ${reffn} \
	${fasta} > ${prefix}.sam    
    samtools view -Sbh ${prefix}.sam > ${prefix}.bam   
    samtools sort -@ ${threads} ${prefix}.bam -o ${prefix}.sorted.bam
    samtools index ${prefix}.sorted.bam
    rm -f ${prefix}.sam
    rm -f ${prefix}.bam
}

function align_with_minimap2_asm20 {
    fasta=$1
    prefix=$2
    reffn=$3
    threads=$4
    minimap2 -x asm20 \
	-t ${threads} --secondary=yes -L -a ${reffn} \
	${fasta} > ${prefix}.sam    
    samtools view -Sbh ${prefix}.sam > ${prefix}.bam   
    samtools sort -@ ${threads} ${prefix}.bam -o ${prefix}.sorted.bam
    samtools index ${prefix}.sorted.bam
    rm -f ${prefix}.sam
    rm -f ${prefix}.bam
}

function align_and_process {
    sample=$1
    dir=$PWD/run_hifiasm
    outdir=${dir}/${sample}/merged_bam/alg_asm20_to_ref_with_secondarySeq
    mkdir $outdir
    minimap2 -x asm20 --secondary=yes -t 12 -L -a ${reffn} ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta > ${outdir}/${sample}.sam
    samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
    samtools sort -@ 10 ${outdir}/${sample}.bam -o ${outdir}/${sample}.sorted.bam
    samtools index ${outdir}/${sample}.sorted.bam
}
function merge_and_rmdup {
    sample=$1
    dir=$PWD/run_hifiasm
    mkdir ${dir}/${sample}/merged_bam
    samtools merge -f ${dir}/${sample}/merged_bam/merged.bam ${dir}/${sample}/break_at_soft_clip/1_asm20_hifi_asm_to_ref.sorted.bam ${dir}/${sample}/break_at_soft_clip/2_asm20_hifi_asm_to_ref.sorted.bam
    samtools sort -@ 12 ${dir}/${sample}/merged_bam/merged.bam -o ${dir}/${sample}/merged_bam/merged.sorted.bam
    samtools index ${dir}/${sample}/merged_bam/merged.sorted.bam
    samtools fasta --reference ${reffn} ${dir}/${sample}/merged_bam/merged.sorted.bam > ${dir}/${sample}/merged_bam/merged_all_reads.fasta
    seqkit rmdup --by-seq ${dir}/${sample}/merged_bam/merged_all_reads.fasta -o ${dir}/${sample}/merged_bam/merged_all_reads.rmdup.fasta
}

#reffn=/home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta
reffn=/home/zmvanw01/git_repos/immune_receptor_genomics/240520/reference.fasta

if [ ! -s ${outdir}/reads.fasta.fai ]
then
    samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
    samtools faidx ${outdir}/reads.fasta
fi

if [ ! -s ${outdir}/hifiasm/asm.bp.hap2.p_ctg.fasta.fai ]
then
    

    if [ ! -f ${outdir}/hifiasm/asm.bp.hap2.p_ctg.gfa ]
    then
		mkdir -p ${outdir}/hifiasm
        hifiasm \
        -o ${outdir}/hifiasm/asm \
        -t ${threads} \
        ${outdir}/reads.fasta
    fi

    for i in p #r
    do
        gfatools \
            gfa2fa \
            ${outdir}/hifiasm/asm.bp.${i}_utg.gfa > \
            ${outdir}/hifiasm/asm.bp.${i}_utg.fasta
    done

    for i in 1 2
    do
        gfatools \
            gfa2fa \
            ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.gfa > \
            ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta

        samtools faidx ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta
    done
fi


for i in 1 2
do
    fn=asm.bp.hap${i}.p_ctg
    if [ ! -s ${outdir}/hifiasm/${fn}_to_ref.sorted.bam.bai ]
    then
	align_with_minimap2_asm20 \
	    ${outdir}/hifiasm/asm.bp.hap${i}.p_ctg.fasta \
	    ${outdir}/hifiasm/${fn}_to_ref \
	    ${reffn} \
	    ${threads}
    fi
done

for i in p #r
do
    fn=asm.bp.${i}_utg
    if [ ! -s ${outdir}/hifiasm/${fn}_to_ref.sorted.bam.bai ]
    then
	align_with_minimap2_asm20 \
	    ${outdir}/hifiasm/${fn}.fasta \
	    ${outdir}/hifiasm/${fn}_to_ref \
	    ${reffn} \
	    ${threads}
    fi
done


mkdir -p ${outdir}/break_at_soft_clip

for i in 1 2
do
    bam=${outdir}/hifiasm/asm.bp.hap${i}.p_ctg_to_ref.sorted.bam

    if [ ! -s ${outdir}/break_at_soft_clip/${i}_hifi_asm_to_ref.sorted.bam ]
    then
        python /home/zmvanw01/git_repos/wasp/hifi-mapping/extract_soft_clip_seq.py \
        ${bam} > ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta

        samtools faidx ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta

        align_with_minimap2 \
        ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta \
        ${outdir}/break_at_soft_clip/${i}_hifi_asm_to_ref \
        ${reffn} \
        ${threads}
    fi

    if [ ! -s ${outdir}/break_at_soft_clip/${i}_asm20_hifi_asm_to_ref.sorted.bam ]
    then
        align_with_minimap2_asm20 \
        ${outdir}/break_at_soft_clip/${i}_hifi_asm.fasta \
        ${outdir}/break_at_soft_clip/${i}_asm20_hifi_asm_to_ref \
        ${reffn} \
        ${threads}
    fi
done
merge_and_rmdup $sample
align_and_process $sample

