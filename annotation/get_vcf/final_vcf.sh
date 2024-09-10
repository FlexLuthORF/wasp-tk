#!/bin/bash
set -e -x

# Define the main directory for VCF files
scratch="$PWD"
outd="$scratch/run_wasp"

function genotype_SV_regions {
    local bam_path="$1"
    local SV_regions_1bp="$2"
    local sample="$3"

    # Ensure the directory exists
    local sample_vcf_dir="${outd}/${sample}/vcfs"
    mkdir -p "$sample_vcf_dir"

    local sample_sv_results="${sample_vcf_dir}/${sample}_sv_genotype_results.txt"

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
    local reffn="$3"
    local num_threads="$4"
    local SV_regions_entire="$5"
    local changeg="$6"
    local anno_config_file="$7"
    local vcfanno="$8"
    local bed_dir="$9"
    
    local sample_vcf_dir="${outd}/${sample}/vcfs"
    mkdir -p "$sample_vcf_dir"

    local of="${sample_vcf_dir}/${sample}"
    bcftools mpileup -B -a QS -Ou -f "${reffn}" -R ${bed_dir}/IG_loci.bed --threads "${num_threads}" "$bam_file" | \
    bcftools call -m -Oz -o "${of}.vcf.gz"

    bcftools index "${of}.vcf.gz"

    if grep -q -P "\t0/1$" "${sample_vcf_dir}/${sample}_sv_genotype_results.txt"; then
        local output_vcf="${sample_vcf_dir}/${sample}_hemi.vcf"
        
        /opt/wasp/conda/bin/python "${changeg}" "${of}.vcf" "${sample_vcf_dir}/${sample}_sv_genotype_results.txt" \
            "${SV_regions_entire}" "${sample}" > "${output_vcf}"

        # Fixing header issues for VCF
        local header_end=$(grep -n '^#CHROM' "$output_vcf" | cut -d ':' -f 1)
        sed -i "1,${header_end}s/^[^#].*//g" "$output_vcf"
        sed -i '1i##fileformat=VCFv4.2' $output_vcf
        sed -i '2i##FILTER=<ID=PASS,Description="All filters passed">' $output_vcf
        sed -i '3i##bcftoolsVersion=1.19+htslib-1.19.1' $output_vcf
        sed -i '/^$/d' "$output_vcf"

        bgzip -c "${output_vcf}" > "${output_vcf}.gz"
        bcftools index "${output_vcf}.gz"
        "${vcfanno}" "${anno_config_file}" "${output_vcf}.gz" > "${sample_vcf_dir}/${sample}_annotated.vcf"
        gzip "${sample_vcf_dir}/${sample}_annotated.vcf"
    else
        "${vcfanno}" "${anno_config_file}" "${of}.vcf.gz" > "${sample_vcf_dir}/${sample}_annotated.vcf"
        gzip "${sample_vcf_dir}/${sample}_annotated.vcf"
    fi
}

# Inputs from the user or script
sample="$1"
bam_file="$2"
reffn="$3"
num_threads="$4"
SV_regions_entire="/opt/wasp/scripts/annotation/KL_SV_regions_entire.bed"
SV_regions_1bp="/opt/wasp/scripts/annotation/SV_regions_1bp.bed"
changeg="/opt/wasp/scripts/annotation/get_vcf/vcf_processing.py"
anno_config_file="/opt/wasp/scripts/annotation/config.toml"
vcfanno="vcfanno"
bed_dir=$5

# Add and index new read group
samtools addreplacerg -r ID:"${sample}" -r SM:"${sample}" -o "${outd}/$sample/${sample}.editRG.bam" "${bam_file}"
samtools index "${outd}/$sample/${sample}.editRG.bam"

# Run functions
genotype_SV_regions "${outd}/$sample/${sample}.editRG.bam" "$SV_regions_1bp" "$sample"
process_vcf "${outd}/$sample/${sample}.editRG.bam" "$sample" "$reffn" "$num_threads" "$SV_regions_entire" "$changeg" "$anno_config_file" "$vcfanno" "$bed_dir"
