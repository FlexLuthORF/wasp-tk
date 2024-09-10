#!/bin/bash

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 path/to/config.cfg sampleID path/to/ccs.bam"
    exit 1
fi

# Source the config file
CONFIG_FILE=$1
if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    echo "Config file not found!"
    exit 1
fi


sample=$2
ccs=$3
container=$4

outdir=$PWD/run_wasp/${sample}
mkdir -p $outdir


if [[ "${ccs}" == *.bam ]]; then
    samtools view ${ccs} | awk '{ print ">"$1"\n"$10 }' > ${outdir}/reads.fasta
else
    cp ${ccs} ${outdir}/reads.fasta
fi
reads="${outdir}/reads.fasta"
# Processing steps
#singularity exec ${container}
bash /opt/wasp/scripts/annotation/create_fofn_from_asm.sh "${outdir}" "${sample}" "${ccs}"
fofn="${outdir}/fofn.tsv"
bash /opt/wasp/scripts/qc/cov.sh "${sample}" "${ccs}" "${reference_fasta}" "${bed_dir}/IG_loci.bed" "${threads}"
bash /opt/wasp/scripts/hifi-mapping/pipeline.sh "${outdir}" "${ccs}" "${threads}" "${sample}" "${reference_fasta}" "${minimap_option}"
/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/process_alleles.py ${sample} ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ${reference_fasta} ${bed_dir} ${allele_ref_dir} ${outdir}
/opt/wasp/conda/bin/python /opt/wasp/scripts/qc/get_asm_stats.py  ${outdir}/merged_bam/merged_all_reads.rmdup.fasta > ${outdir}/merged_bam/${sample}.asm.stats
samtools stats ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.stats
samtools flagstat ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam > ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.asm-to-ref.flagstats
bash /opt/wasp/scripts/annotation/read-support/get_read_support_VDJs.sh ${fofn} ${reference_fasta} ${bed_dir}/IG_loci.bed ${threads} ${outdir} ${ccs_minimap_option}
bash /opt/wasp/scripts/annotation/get_vcf/final_vcf.sh ${sample} ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ${reference_fasta} ${threads} ${bed_dir}
