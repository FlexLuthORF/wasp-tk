#!/bin/bash
set -e -x

outdir=$1
ccs=$2
threads=$3

function align_with_minimap2 {
    fasta=$1
    prefix=$2
    reffn=$3
    threads=$4
    minimap2 \
	-t ${threads} -L -a ${reffn} \
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
	-t ${threads} -L -a ${reffn} \
	${fasta} > ${prefix}.sam    
    samtools view -Sbh ${prefix}.sam > ${prefix}.bam   
    samtools sort -@ ${threads} ${prefix}.bam -o ${prefix}.sorted.bam
    samtools index ${prefix}.sorted.bam
    rm -f ${prefix}.sam
    rm -f ${prefix}.bam
}

reffn=../../immune_receptor_genomics/current/reference.fasta

if [ ! -s ${outdir}/reads.fasta.fai ]
then
    samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
    samtools faidx ${outdir}/reads.fasta
fi

if [ ! -s ${outdir}/hifiasm/asm.bp.hap2.p_ctg.fasta.fai ]
then
    mkdir -p ${outdir}/hifiasm
    hifiasm \
	-o ${outdir}/hifiasm/asm \
	-t ${threads} \
	${outdir}/reads.fasta
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


for iter in 1 2
do
    mkdir -p ${outdir}/break_at_soft_clip/${iter}
    for i in 1 2
    do
	if [ ${iter} == 1 ]
	then
	    bam=${outdir}/hifiasm/asm.bp.hap${i}.p_ctg_to_ref.sorted.bam
	else
	    bam=${outdir}/break_at_soft_clip/1/${i}_hifi_asm_to_ref.sorted.bam
	fi
	if  [ ! -s ${outdir}/break_at_soft_clip/${iter}/${i}_hifi_asm_to_ref.sorted.bam ]
	then
	    python extract_soft_clip_seq.py \
		${bam} > ${outdir}/break_at_soft_clip/${iter}/${i}_hifi_asm.fasta
	    samtools faidx ${outdir}/break_at_soft_clip/${iter}/${i}_hifi_asm.fasta
	    align_with_minimap2 \>
		${outdir}/break_at_soft_clip/${iter}/${i}_hifi_asm.fasta \
		${outdir}/break_at_soft_clip/${iter}/${i}_hifi_asm_to_ref \
		${reffn} \
		${threads}
	fi
	if [ ! -s ${outdir}/break_at_soft_clip/${iter}/${i}_asm20_hifi_asm_to_ref.sorted.bam ]
	then
	    align_with_minimap2_asm20 \
		${outdir}/break_at_soft_clip/${iter}/${i}_hifi_asm.fasta \
		${outdir}/break_at_soft_clip/${iter}/${i}_asm20_hifi_asm_to_ref \
		${reffn} \
		${threads}
	fi
    done
done
