#!/bin/bash

outdir="$1"
sample="$2"
ccs="$3"
fofn="${outdir}/fofn.tsv"
asm_bam="${outdir}/merged_bam/final_asm20_to_ref_with_secondarySeq/${sample}.sorted.bam"

# Initialize the line with sample, asm_bam, and ccs_bam
line="${sample}\t${asm_bam}"

# Define locus mappings
declare -A loci=(
    [chr2_gene]='IGK'
    [chr2_import]='IGK'
    [chr22_gene]='IGL'
    [chr22_import]='IGL'
    [igh_gene]='IGH'
    [igh_import]='IGH'
    [ighc_gene]='IGHC'
    [ighc_import]='IGHC'
    [trb_gene]='TRB'
    [trb_import]='TRB'
    [trg_gene]='TRG'
    [trg_import]='TRG'
    [trd_gene]='TRD'
    [trd_import]='TRD'
    [tra_gene]='TRA'
    [tra_import]='TRA'
)

# Loop through loci and append gene and import files
for key in "${!loci[@]}"; do
    dir="${loci[$key]}"
    gene_file="${outdir}/merged_bam/annotations/${sample}/${dir}/${sample}__make_gene_file.csv"
    import_file="${outdir}/merged_bam/annotations/${sample}/${dir}/${sample}_make_gene_file_imported.csv"
    
    line+="\t${gene_file}\t${import_file}"
done

# Append ccs_bam at the end
line+="\t${ccs}"

# Output to fofn.tsv
echo -e "$line" > "$fofn"
