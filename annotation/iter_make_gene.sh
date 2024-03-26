#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 path/to/your/fofn.tsv"
    exit 1
fi

fofn_path="$1"
base_output_dir="$PWD"

# Correctly mapping loci to their chromosome names for the command, but directories are named igh, igk, igl
declare -A loci=( ["IGH"]="igh" ["IGK"]="igk" ["IGL"]="igl" )
declare -A chroms=( ["IGH"]="igh" ["IGK"]="chr2" ["IGL"]="chr22" )

for locus in "${!loci[@]}"; do
    strand="+"
    if [ "$locus" == "IGL" ]; then
        strand="+"
    else
        strand="\"+-\""
    fi

    # Ensure the specific output directory exists
    mkdir -p "${base_output_dir}/${loci[$locus]}"

    command="python /home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/make_gene_file.py -l ${locus} -r ${chroms[$locus]} -s ${strand} -b /home/zmvanw01/git_repos/immune_receptor_genomics/current -f /home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta -o ${base_output_dir}/${loci[$locus]} -n ${fofn_path}"
    sbatch_cmd="sbatch --wrap=\"${command}\""
    eval $sbatch_cmd
done

echo "Submission completed."
