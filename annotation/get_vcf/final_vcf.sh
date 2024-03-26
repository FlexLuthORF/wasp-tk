#!/bin/bash
set -e -x

# Function to genotype SV regions
function genotype_SV_regions {
    local bam_path="$1"
    local SV_regions_1bp="$2"
    local outd="$3"
    local sample="$4"

    sample_sv_results="${outd}/${sample}_sv_genotype_results.txt"

    while read -r sv_region; do
        sv_name=$(echo "$sv_region" | cut -f4)
        grep -w "$sv_name" "$SV_regions_1bp" > "${sv_name}.bed"

        mpileup_output=$(samtools mpileup -l "${sv_name}.bed" -f "${reffn}" "$bam_path" | head -1)
        genotype=$(echo "$mpileup_output" | awk '{print $5}')
        asterisk_count=$(echo "$genotype" | tr -cd '*' | wc -c)
        total_count=$(echo "$genotype" | wc -c)

        if [ "$asterisk_count" -eq 0 ]; then
            genotype_label="0/0"
        elif [ "$asterisk_count" -eq "$((total_count-1))" ]; then
            genotype_label="1/1"
        else
            genotype_label="0/1"
        fi

        echo -e "$sample\t$sv_name\t$genotype_label" >> "$sample_sv_results"
    done < "$SV_regions_1bp"
}

function process_vcf {
    local bam_file="$1"
    local sample="$2"
    local outd="$3"
    local reffn="$4"
    local num_threads="$5"
    local SV_regions_entire="$6"
    local changeg="$7"
    local anno_config_file="$8"
    local vcfanno="$9"

    sample_outd="${outd}/${sample}"
    mkdir -p "${sample_outd}"

    of="${sample_outd}/${sample}"
    bcftools mpileup -B -a QS -Ou -f "${reffn}" \
        --threads "${num_threads}" "$bam_file" | \
        bcftools call -m -Oz -o "${of}.vcf.gz"
    bcftools index "${of}.vcf.gz"

    bcftools mpileup -B -a QS -Ou -f "${reffn}" \
        --threads "${num_threads}" "$bam_file" | \
        bcftools call -m -Ov -o "${of}.vcf"

    mkdir -p "${sample_outd}/change_to_hemi"
    mkdir -p "${outd}/annotated_vcfs/${sample}"

    sample_sv_results="${outd}/${sample}_sv_genotype_results.txt"

    # Check if the sample needs hemizygous adjustment
    if grep -q -P "\t0/1$" "$sample_sv_results"; then
        # The sample has at least one 0/1 genotype, so it needs adjustment
        output_vcf="${sample_outd}/change_to_hemi/${sample}_hemi.vcf"
        python "${changeg}" "${of}.vcf" "$sample_sv_results" \
            "${SV_regions_entire}" "${sample}" > "${output_vcf}"
        
    # THIS FIXES A WEIRD CYVCF2 ISSUE WE ARE HAVING
    # Remove lines at the start of the file that don't start with "##"
        header_end=$(grep -n '^#CHROM' "$output_vcf" | cut -d ':' -f 1)
        sed -i "1,${header_end}s/^[^#].*//g" "$output_vcf"
        # Add the required header lines at the beginning of the file
        sed -i '1i##fileformat=VCFv4.2' $output_vcf
        sed -i '2i##FILTER=<ID=PASS,Description="All filters passed">' $output_vcf
        sed -i '3i##bcftoolsVersion=1.19+htslib-1.19.1' $output_vcf
        # Remove empty lines from the file
        sed -i '/^$/d' "$output_vcf"

        bgzip -c "${output_vcf}" > "${output_vcf}.gz"
        bcftools index "${output_vcf}.gz"
        "${vcfanno}" "${anno_config_file}" "${output_vcf}.gz" \
            > "${outd}/annotated_vcfs/${sample}/${sample}_annotated.vcf"
    else
        # The sample doesn't need hemizygous adjustment
        "${vcfanno}" "${anno_config_file}" "${of}.vcf.gz" > "${outd}/annotated_vcfs/${sample}/${sample}_annotated.vcf"
    fi
}

# Main script
scratch="$PWD"
outd="${scratch}/geno_analysis/per_samp"
mkdir -p "${outd}"

bam_file="$1"
sample=$(basename "$bam_file" .bam | cut -d '_' -f 1)
reffn="$2"
num_threads="$3"
SV_regions_entire="$4"
SV_regions_1bp="$5"
changeg="/home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/get_vcf/vcf_processing.py"
anno_config_file="$6"
vcfanno="$7"

samtools addreplacerg -r ID:"${sample}" -r SM:"${sample}" \
    -o "${scratch}/$sample/${sample}.editRG.bam" "${bam_file}"
samtools index "${scratch}/$sample/${sample}.editRG.bam"

bam_path="${scratch}/$sample/${sample}.editRG.bam"

genotype_SV_regions "$bam_path" "$SV_regions_1bp" "$outd" "$sample"
process_vcf "$bam_path" "$sample" "$outd" "$reffn" "$num_threads" "$SV_regions_entire" "$changeg" "$anno_config_file" "$vcfanno"