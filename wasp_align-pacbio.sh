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

outdir=$PWD/run_hifiasm/${sample}
mkdir -p $outdir

# Processing steps
#singularity exec ${container}
bash /opt/wasp/scripts/annotation/create_fofn_from_asm.sh "${outdir}" "${sample}" "${ccs}"
fofn="${outdir}/fofn.tsv"
bash /opt/wasp/scripts/hifi-mapping/pipeline.sh ${outdir} ${ccs} ${threads} ${sample} ${reference_fasta}
/opt/wasp/conda/binpython /opt/wasp/scripts/annotation/process_alleles.py ${sample} ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ${reference_fasta} ${bed_dir} ${allele_ref_dir} ${outdir}
bash /opt/wasp/scripts/annotation/read-support/get_read_support_VDJs.sh ${fofn} ${reference_fasta} ${bed_dir}/IG_loci.bed ${threads}
