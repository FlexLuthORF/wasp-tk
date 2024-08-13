#!/bin/bash
set -e -x

scratch=$PWD
sample=$1
ccs=$2
reffn=$3
#refbed=/home/zmvanw01/test-beds/sorted_region.bed
#refbed=/home/zmvanw01/projects/t_Parks/parks.bed
#refbed=/home/zmvanw01/240520-coverage.bed
refbed=$4
threads=$5

function run_get_ccs_cov {
    mkdir -p ${scratch}/run_wasp/${sample}/ccs_cov
    outd=${scratch}/run_wasp/${sample}/ccs_cov
    bam_path=${scratch}/run_wasp/${sample}/ccs_cov/ccs_to_ref.sorted.bam
    pers_bam_path=${scratch}/run_wasp/${sample}/ccs_cov/ccs_to_pers-ref.sorted.bam
    #bamtocounts --coords ${refbed} ${bam_path} > ${outd}/${sample}/${sample}_ccs-cov_counts.bed
    bamtocov --regions ${refbed} --report ${outd}/$sample/${sample}_stats.tsv ${bam_path} > ${outd}/$sample/${sample}_cov-cov.bed
    #bamtocov --regions ${refbed} --report ${outd}/$sample/${sample}_stats_pers-ref.tsv ${pers_bam_path} > ${outd}/$sample/${sample}_cov-cov_pers-ref.bed
}
function map_ccs_to_ref {
    dir=$scratch/ccs_cov
    outdir=${scratch}/run_wasp/${sample}/ccs_cov
    mkdir -p $outdir
    samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
    samtools faidx ${outdir}/reads.fasta
    minimap2 -x map-hifi --secondary=no -t "${threads}" -L -a ${reffn} ${outdir}/reads.fasta > ${outdir}/${sample}.sam
    samtools view -Sbh ${outdir}/${sample}.sam > ${outdir}/${sample}.bam
    samtools sort -@ 12 ${outdir}/${sample}.bam -o ${outdir}/ccs_to_ref.sorted.bam
    samtools index ${outdir}/ccs_to_ref.sorted.bam
}
function map_ccs_to_pers_ref {
    dir=$scratch/ccs_cov
    outdir=${scratch}/run_wasp/${sample}/ccs_cov
    mkdir -p $outdir
    minimap2 -x map-hifi --secondary=no -t 12 -L -a ${reffn} $PWD/run_wasp/${sample}/read_support/${sample}/ccs_to_pers/pers_ref.fasta > ${outdir}/${sample}.pers.sam
    samtools view -Sbh ${outdir}/${sample}.pers.sam > ${outdir}/${sample}.pers.bam
    samtools sort -@ 12 ${outdir}/${sample}.pers.bam -o ${outdir}/ccs_to_pers-ref.sorted.bam
    samtools index ${outdir}/ccs_to_pers-ref.sorted.bam
}
#map_ccs_to_pers_ref
map_ccs_to_ref
run_get_ccs_cov
