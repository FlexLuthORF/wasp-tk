#!/bin/bash

outdir="$1"
sample="$2"
ccs="$3"
asm_bam="${4:-${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam}"
fofn="${outdir}/fofn.tsv"

# Initialize the line with sample, asm_bam
line="${sample}\t${asm_bam}"

# Define locus mappings and their order
declare -a order=(chr2_gene chr22_gene igh_gene ighc_gene trb_gene trg_gene trd_gene tra_gene)
declare -A loci=(
    [chr2_gene]='IGK'
    [chr22_gene]='IGL'
    [igh_gene]='IGH'
    [ighc_gene]='IGHC'
    [trb_gene]='TRB'
    [trg_gene]='TRG'
    [trd_gene]='TRD'
    [tra_gene]='TRA'
)

# Loop through loci in specified order and append gene and import files
for key in "${order[@]}"; do
    dir="${loci[$key]}"
    gene_file="${outdir}/annotations/${sample}/${dir}/${sample}_make_gene_file.csv"
    import_file="${outdir}/annotations/${sample}/${dir}/${sample}_make_gene_file_imported.csv"
    
    line+="\t${gene_file}\t${import_file}"
done

# Append ccs_bam at the end
line+="\t${ccs}"

# Output to fofn.tsv
echo -e "$line" > "$fofn"
