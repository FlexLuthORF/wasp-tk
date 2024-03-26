#!/bin/bash
# Usage:
# ./script_name.sh scratch changeg vcfanno reffn fofn.tsv num_threads SV_regions_entire SV_regions_1bp anno_config_file IG_loci

set -e -x

# Assigning command line arguments to variables
scratch="$1" # can be $pwd
changeg="$2"
vcfanno="$3"
reffn="$4"
fofn="$5"
num_threads="$6"
SV_regions_entire="$7"
SV_regions_1bp="$8"
anno_config_file="$9"

# Function for genotyping the SV regions
function genotype_SV_regions_entire {
    mkdir -p "${scratch}/geno_analysis/per_samp"
    outd="${scratch}/geno_analysis/per_samp"
    > "${outd}/SV_genotype_results.txt"

    while IFS=$'\t' read -r sample bam_file_path; do
        

        samtools addreplacerg -r ID:"${sample}" -r SM:"${sample}" \
            -o "${scratch}/$sample/${bam_base_name}.editRG.bam" "${bam_file_path}"
        samtools index "${scratch}/$sample/${bam_base_name}.editRG.bam"

        bam_base_name=$(basename "${bam_file_path}" .bam)
        bam_path="${scratch}/$sample/${bam_base_name}.editRG.bam"
        echo $bam_path

        grep "IGKV1-NL1" $SV_regions_1bp > IGKV1NL1.bed 
        mpileup_output_IGKV1NL1=$(samtools mpileup -l IGKV1NL1.bed "$bam_path" | head -1)
        genotype_IGKV1NL1=$(echo "$mpileup_output_IGKV1NL1" | awk '{print $5}')
        asterisk_count_IGKV1NL1=$(echo "$genotype_IGKV1NL1" | tr -cd '*' | wc -c)
        
        if [ "$asterisk_count_IGKV1NL1" -eq 0 ]; then
            genotype_label_IGKV1NL1="0/0"
        elif [ "$asterisk_count_IGKV1NL1" -eq 1 ]; then
            genotype_label_IGKV1NL1="0/1"
        elif [ "$asterisk_count_IGKV1NL1" -eq 2 ]; then
            genotype_label_IGKV1NL1="1/1"
        else
            genotype_label_IGKV1NL1="unknown"
        fi
        
        grep "IGLV5-39" $SV_regions_1bp > IGLV539.bed
        mpileup_output_IGLV539=$(samtools mpileup -l IGLV539.bed -f "${reffn}" "$bam_path"| head -1)
        genotype_IGLV539=$(echo "$mpileup_output_IGLV539" | awk '{print $5}')
        asterisk_count_IGLV539=$(echo "$genotype_IGLV539" | tr -cd '*' | wc -c)
        
        if [ "$asterisk_count_IGLV539" -eq 0 ]; then
            genotype_label_IGLV539="0/0"
        elif [ "$asterisk_count_IGLV539" -eq 1 ]; then
            genotype_label_IGLV539="0/1"
        elif [ "$asterisk_count_IGLV539" -eq 2 ]; then
            genotype_label_IGLV539="1/1"
        else
            genotype_label_IGLV539="unknown"
        fi

        echo -e "$sample\t$genotype_label_IGKV1NL1\t$genotype_label_IGLV539" >> "${outd}/SV_genotype_results.txt"
    done < "${fofn}"
}

function get_vcf_files {
    mkdir -p "${scratch}/geno_analysis/per_samp"
    outd="${scratch}/geno_analysis/per_samp"

    cat ${fofn} | while read sample bam_file_path; do
        echo "Reading sample: $sample, BAM file path: $bam_file_path"  # Debug print

        bam_dir=$(dirname "${bam_file_path}")
        echo "BAM directory: $bam_dir"  # Debug print

        bam_file="${bam_dir}/${sample}_merged.sorted.editRG.bam"
        echo "BAM file to be processed: $bam_file"  # Debug print

        mkdir -p "${outd}/${sample}"
        echo "Output directory created for sample: ${outd}/${sample}"  # Debug print

        of="${outd}/${sample}/${sample}"
        echo "Output file path: $of"  # Debug print

       
        #here
        
        bcftools mpileup -B -a QS -Ou -f "${reffn}" \
            --threads "${num_threads}" "$bam_file" | \
            bcftools call -m -Oz -o "${of}.vcf.gz"
        bcftools index "${of}.vcf.gz"

        bcftools mpileup -B -a QS -Ou -f "${reffn}" \
            --threads "${num_threads}" "$bam_file" | \
            bcftools call -m -Ov -o "${of}.vcf"
        

        
        #python "${changeg}" "${of}.vcf" "${SV_regions_entire}" "${outd}/${sample}/change_to_hemi/${sample}.vcf"
       # bgzip -c "${outd}/${sample}/change_to_hemi/${sample}.vcf" > "${outd}/${sample}/change_to_hemi/${sample}.vcf.gz"
        #bcftools index "${outd}/${sample}/change_to_hemi/${sample}.vcf.gz"
        #
        #"${vcfanno}" "${anno_config_file}" "${outd}/${sample}/change_to_hemi/${sample}.vcf.gz" > "${outd}/annotated_vcfs/${sample}/${sample}_annotated.vcf"
    done

    cat "${outd}/SV_genotype_results.txt" | while read sample IGK_SV_GT IGL_SV_GT; do   
        of="${outd}/${sample}/${sample}"
        mkdir -p "${outd}/${sample}/change_to_hemi"
        mkdir -p "${outd}/annotated_vcfs/${sample}"
        
        if [ "$IGK_SV_GT" == "0/1" ] && [ "$IGL_SV_GT" == "0/1" ]; then
        	sv_regions_input="${SV_regions_entire}"
    	elif [ "$IGK_SV_GT" == "0/1" ] && [ "$IGL_SV_GT" != "0/1" ]; then
        	grep "IGKV1-NL1" "${SV_regions_entire}" > IGKV1-NL1.bed
            sv_regions_input=IGKV1-NL1.bed  
    	elif [ "$IGK_SV_GT" != "0/1" ] && [ "$IGL_SV_GT" == "0/1" ]; then
        	grep "IGLV5-39" "${SV_regions_entire}" > IGLV5-39.bed
            sv_regions_input=IGLV5-39.bed
    	else
        	continue # Skip to the next iteration without executing the python script
    	fi
        python "${changeg}" "${of}.vcf" \
            "${sv_regions_input}" \
            "${outd}/${sample}/change_to_hemi/${sample}.vcf"
        bgzip -c "${outd}/${sample}/change_to_hemi/${sample}.vcf" > "${outd}/${sample}/change_to_hemi/${sample}.vcf.gz"
        bcftools index "${outd}/${sample}/change_to_hemi/${sample}.vcf.gz"
        mkdir -p "${outd}/annotated_vcfs/${sample}"
        "${vcfanno}" "${anno_config_file}" "${outd}/${sample}/change_to_hemi/${sample}.vcf.gz" \
            > "${outd}/annotated_vcfs/${sample}/${sample}_annotated.vcf"
    done 


#annotate the rest of them
    cat "${outd}/SV_genotype_results.txt" | while read sample IGK_SV_GT IGL_SV_GT; do
        if [ "$IGK_SV_GT" != "0/1" ] || [ "$IGL_SV_GT" != "0/1" ]; then
            of="${outd}/${sample}/${sample}"
            mkdir -p "${outd}/annotated_vcfs/${sample}"
            "${vcfanno}" "${anno_config_file}" "${of}.vcf.gz" > "${outd}/annotated_vcfs/${sample}/${sample}_annotated.vcf"
        fi
    done
}

genotype_SV_regions_entire
get_vcf_files
