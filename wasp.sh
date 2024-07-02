#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 path/to/config.cfg path/to/input_fofn.tsv"
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

user=$(whoami)
input_fofn=$2

cat $input_fofn | while read sample ccs
do
    outdir=$PWD/run_hifiasm/${sample}
    mkdir -p $outdir
    # Create the sbatch script for each sample
    echo "#!/bin/bash
#SBATCH --time=88:00:00
#SBATCH -p compute
#SBATCH -o ${outdir}/job.txt

mkdir -p $outdir
bash create_fofn.sh '${outdir}' '${sample}' '${ccs}'
fofn='${outdir}/fofn.tsv'
bash pipeline.sh ${outdir} ${ccs} ${threads} ${sample} ${reference_fasta}
python process_alleles.py ${sample} ${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam ${reference_fasta} ${bed_dir} ${allele_ref_dir}
bash get_read_support_VDJs.sh ${fofn} ${reference_fasta} ${bed_dir}/IG_loci.bed" > ${outdir}/${sample}.sbatch

    # Submit the job
    sbatch ${outdir}/${sample}.sbatch
    echo "${sample} job submitted"


    # Check queue and pause if there are too many jobs running
    count=$(squeue | grep $user | wc -l)
    while [ $count -gt $((max_jobs - 1)) ]
    do
        sleep 20s
        count=$(squeue | grep $user | wc -l)
    done

done

