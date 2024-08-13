#!/bin/bash

# Check arguments
#if [ "$#" -ne 4 ]; then
 #   echo "Usage: $0 path/to/config.cfg sampleID path/to/ccs.bam path/to/asm.bam"
  #  exit 1
#fi

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
asm=$4
outdir=${5:-$PWD/run_wasp/${sample}}

mkdir -p $outdir

# Processing steps
#singularity exec ${container}
bash /opt/wasp/scripts/annotation/create_fofn_from_asm.sh "${outdir}" "${sample}" "${ccs}" "${asm}"
fofn="${outdir}/fofn.tsv"
/opt/wasp/conda/bin/python /opt/wasp/scripts/annotation/process_alleles.py ${sample} ${asm} ${reference_fasta} ${bed_dir} ${allele_ref_dir} ${outdir}
bash /opt/wasp/scripts/annotation/read-support/get_read_support_VDJs.sh ${fofn} ${reference_fasta} ${bed_dir}/IG_loci.bed ${threads} ${outdir} ${ccs_minimap_option}
