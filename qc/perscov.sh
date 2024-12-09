#!/bin/bash

# Check if correct number of arguments was given
# if [ "$#" -ne 3 ]; then
#     echo "Usage: $0 <ccs_reads_bam> <contigs_to_ref_bam> <IG_loci.bed>"
#     exit 1
# fi

# Assign input arguments to variables
sample="$1"
ccs_reads_bam="$2"
contigs_to_ref_bam="$3"
ig_loci_bed="$4"
outdir="$5"

# Step 1: Calculate coverage for each contig using bamtocov
bamtocov --wig 0 -o ${outdir}/ccs_cov/coverage_output.bed "$ccs_reads_bam"

# Step 2: Extract contig mappings from the other BAM and convert to BED format
samtools view "$contigs_to_ref_bam" | awk 'BEGIN {OFS="\t"} {print $3, $4, $4+length($10), $1}' > ${outdir}/ccs_cov/contig_to_chrom_coords.bed

# Step 3: Filter contigs that overlap with IG loci using bedtools
bedtools intersect -a ${outdir}/ccs_cov/contig_to_chrom_coords.bed -b "$ig_loci_bed" -wo > ${outdir}/ccs_cov/filtered_contigs.bed

# Step 4: Combine coverage with mapping information from filtered contigs and calculate average coverage per chromosome
awk 'NR==FNR{cov[$1]=$3; next} {print $1, cov[$4]}' ${outdir}/ccs_cov/coverage_output.bed ${outdir}/ccs_cov/filtered_contigs.bed | \
awk 'BEGIN {OFS="\t"} {if($2!="") {chr[$1]+=$2; count[$1]++}} END {for (c in chr) print c, chr[c]/count[c]}' > ${outdir}/ccs_cov/average_chrom_coverage.tsv

echo "Coverage calculation and mapping complete."