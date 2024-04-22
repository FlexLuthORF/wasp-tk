#!/bin/bash
set -e -x

reffn=/home/egenge01/projects/IGL_ref_mod/reference_ready/modified_reference_renamed.fasta
IG_loci=/home/egenge01/projects/12_sample_test/IG_loci.bed
masked_ref=/home/egenge01/projects/12_sample_test/reference_IGloci_masked.fasta
scratch=$PWD
mask_ref=${scratch}/ref_IG_masked.fasta

function run_make_ref_masked {

bedtools maskfasta -fi ${reffn} -bed ${IG_loci} -fo ${scratch}/ref_IG_masked.fasta
samtools faidx ${scratch}/ref_IG_masked.fasta

}

function run_map_ccs_to_pers {
    cat extended_samples_paths.fofn | while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import ccs_bam
    do
	mkdir -p ${scratch}/read_support/${sample}
	outd=${scratch}/read_support/${sample}
	mkdir -p ${outd}/ccs_to_pers
	#convert pacbio hifi reads to fasta
	samtools view ${ccs_bam} | awk '{ print ">"$1"\n"$10 }' > ${outd}/ccs_to_pers/reads.fasta
	samtools faidx ${outd}/ccs_to_pers/reads.fasta
	
	#create personalized reference
	samtools view -F 0x100 -F 0x800 "${asm_bam}" | awk '{print ">"$1"\n"$10}' > "${outd}/ccs_to_pers/contigs.fasta"
	samtools faidx ${outd}/ccs_to_pers/contigs.fasta
	pers_contigs=${outd}/ccs_to_pers/contigs.fasta
	cat ${mask_ref} ${pers_contigs} \
	    > ${outd}/ccs_to_pers/pers_ref.fasta
	samtools faidx ${outd}/ccs_to_pers/pers_ref.fasta
	wait
	sbatch --time=88:00:00 -p compute -o ${outd}/ccs_to_pers/job_map_to_pers.txt --wrap="/home/egenge01/minimap2/minimap2 -ax map-hifi --secondary=no -t 10 -L ${outd}/ccs_to_pers/pers_ref.fasta ${outd}/ccs_to_pers/reads.fasta > ${outd}/ccs_to_pers/output.sam;
 samtools view -Sbh ${outd}/ccs_to_pers/output.sam > ${outd}/ccs_to_pers/output.bam;
 samtools sort -@ 10 ${outd}/ccs_to_pers/output.bam -o ${outd}/ccs_to_pers/output.sorted.bam;
 samtools index ${outd}/ccs_to_pers/output.sorted.bam;
 wait
  rm -f ${outd}/ccs_to_pers/output.sam"
    done
}

function run_append_pos {
    while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import ccs_bam
    do
        # Base directory for imported genes
        base_outd="${scratch}/read_support/${sample}/imported_genes"
        
        # Create subdirectories for each type and define output paths
        mkdir -p "${base_outd}/chr2"
        chr2_import_out="${base_outd}/chr2/$(basename "${chr2_import}")"

        mkdir -p "${base_outd}/chr22"
        chr22_import_out="${base_outd}/chr22/$(basename "${chr22_import}")"

        mkdir -p "${base_outd}/igh"
        igh_import_out="${base_outd}/igh/$(basename "${igh_import}")"

#        mkdir -p "${base_outd}/ighc"
#        ighc_import_out="${base_outd}/ighc/$(basename "${ighc_import}")"

        # Call the Python script with input, original output, and modified output paths
        python append_pos_import_genes.py "${chr2_gene}" "${chr2_import}" "${chr2_import_out}"
        python append_pos_import_genes.py "${chr22_gene}" "${chr22_import}" "${chr22_import_out}"
        python append_pos_import_genes.py "${igh_gene}" "${igh_import}" "${igh_import_out}"
##      python append_pos_import_genes.py "${ighc_gene}" "${ighc_import}" "${ighc_import_out}"
#        python append_pos_import_genes_ighc.py "${ighc_gene}" "${ighc_import}" "${ighc_import_out}"
    done < extended_samples_paths.fofn
}

function get_read_support_vdj3 {
    while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import ccs_bam
    do
        base_outd="${scratch}/read_support/${sample}/imported_genes"
        bam_file="${scratch}/read_support/${sample}/ccs_to_pers/output.sorted.bam" # ccs reads to personalized reference
        ref="${scratch}/read_support/${sample}/ccs_to_pers/pers_ref.fasta" # personalized reference

        if [ ! -f "${bam_file}.bai" ]; then
            samtools index "$bam_file"
        fi

        for gene_type in "chr2" "chr22" "igh" #"ighc"
        do
            import_out="${base_outd}/${gene_type}/${sample}_make_gene_file_imported.csv"

            if [[ -f "$import_out" ]]; then
		modified_import_out="$import_out"
		
                tmp_file="${import_out}_read_support.tmp"
                echo "Total_Positions,Average_Coverage,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Percent_Accuracy,Positions_With_At_Least_10x_Coverage,Fully_Spanning_Reads,Fully_Spanning_Reads_100%_Match" > "$tmp_file"
                header=$(head -n 1 "$import_out")
                IFS=',' read -ra header_cols <<< "$header"
                for i in "${!header_cols[@]}"; do
                    case "${header_cols[$i]}" in
                        "contig") contig_col=$i ;;
                        "REGION_start") start_col=$i ;;
                        "REGION_end") end_col=$i ;;
			"gene") gene_col=$i ;;  # Added case for 'gene'
                    esac
                done
		tmp_counts="${tmp_file}_counts"  # Temporary file to hold counts
		> "$tmp_counts"  # Clear or create the temp file for counts

                tail -n +2 "$import_out" | while IFS=, read -ra line
                do
                    contig="${line[$contig_col]}"
                    start=$(echo "${line[$start_col]}" | awk '{printf "%.0f", $1}')
                    end=$(echo "${line[$end_col]}" | awk '{printf "%.0f", $1}')
                    
		    gene="${line[$gene_col]}"  # Extract the 'gene' value using the identified column index
                    region="${contig}:${start}-${end}"

                    contig_filename=$(echo "$contig" | tr '/' '_')
                    tmp_bam="${base_outd}/${gene_type}/${contig_filename}_${start}_${end}.bam"
                    mkdir -p "$(dirname "$tmp_bam")"
                    samtools view -F 0x100 -F 0x800 -b "$bam_file" -o "$tmp_bam" -U "/dev/null" "${contig}:${start}-${end}"
                    samtools index "$tmp_bam"

		    samtools mpileup -f "$ref" -r "$region" "$tmp_bam" | \
			awk -v total_positions="$((end - start + 1))" -v sample="$sample" \
			'BEGIN {
    total_reads=0; mismatched_positions=0; matched_positions=0; positions_with_10x=0;
    mismatch_list=""; match_list="";
}
{
    total_reads += length($5);
    mismatches = length(gensub(/[.,]/, "", "g", $5));
    matches = length(gensub(/[^.,]/, "", "g", $5));
    mismatch_list = (mismatch_list == "" ? mismatches : mismatch_list ":" mismatches);
    match_list = (match_list == "" ? matches : match_list ":" matches);

    coverage = length($5);
    if (coverage >= 10) {
        positions_with_10x++;
    }

    mismatch_rate = mismatches / coverage;
    match_rate = matches / coverage;

    if (mismatch_rate > 0.2) {
        mismatched_positions++;
    }

    if (match_rate > 0.8) {
        matched_positions++;
    }
}
END {
    avg_reads_per_position = (total_positions > 0) ? total_reads / total_positions : 0;
    percent_accuracy = (matched_positions / total_positions) * 100;
    print total_positions, avg_reads_per_position, mismatched_positions, matched_positions, mismatch_list, match_list, percent_accuracy, positions_with_10x;}' OFS=',' >> "${tmp_file}_awk_out"
		    python match_subsequences3.py "$tmp_bam" "$contig" "$start" "$end" "$gene" "$import_out" > "${tmp_file}_py_out"
		    wait
		    paste -d ',' "${tmp_file}_awk_out" "${tmp_file}_py_out" >> "$tmp_file"
		    rm "${tmp_file}_awk_out" "${tmp_file}_py_out"
                    rm "$tmp_bam" "${tmp_bam}.bai"
                done

# New block to merge data and update Subject and Sample_Name
		if [[ -f "$import_out" && -f "$tmp_file" ]]; then
		    combined_file="${import_out%.csv}_combined.csv"
		    final_output="${import_out%.csv}_with_read_support.csv"
		    
		    # Merge the original import file with the new tmp file
		    paste -d ',' "$import_out" "$tmp_file" > "$combined_file"
		    rm -f "$tmp_file"
		    
		    # Update the Subject and Sample_Name columns in the final output
		    awk -v sample="$sample" 'BEGIN{FS=OFS=","} {
    if (NR == 1) {
        print;
    } else {
        $2 = sample;  # Update the 2nd column with $sample
        $3 = sample;  # Update the 3rd column with $sample
        print;
    }
}' "$combined_file" > "$final_output"

		    rm -f "$combined_file"
		fi
	    fi
        done
    done < extended_samples_paths.fofn
}


function get_read_support_ighc {
    while read sample asm_bam chr2_gene chr2_import chr22_gene chr22_import igh_gene igh_import ighc_gene ighc_import ccs_bam
    do
        base_outd="${scratch}/read_support/${sample}/imported_genes"
        bam_file="${scratch}/read_support/${sample}/ccs_to_pers/output.sorted.bam" # ccs reads to personalized reference
        ref="${scratch}/read_support/${sample}/ccs_to_pers/pers_ref.fasta" # personalized reference

        if [ ! -f "${bam_file}.bai" ]; then
            samtools index "$bam_file"
        fi

        for gene_type in "ighc" #"chr2" "chr22" "igh"
        do
            import_out="${base_outd}/${gene_type}/${sample}_make_gene_file_imported.csv"

            if [[ -f "$import_out" ]]; then
                # Use csvcut to remove the "notes" column, which causes problems and can have or always has a comma
                csvcut -C "notes" "$import_out" > "${import_out%.csv}_nonotes.csv"
                modified_import_out="${import_out%.csv}_nonotes.csv"

                tmp_file="${import_out}_read_support.tmp"
                echo "Total_Positions,Total_Reads_by_Positions,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Mismatched_Positions_Coverage_Less_Than_10,Mismatched_Positions_Coverage_10_Or_Greater,Matched_Positions_Coverage_Less_Than_10,Matched_Positions_Coverage_10_Or_Greater,Percent_Accuracy" > "$tmp_file"
		
		header=$(head -n 1 "$modified_import_out")
		IFS=',' read -ra header_cols <<< "$header"
		declare -A exon_cols  # Use an associative array to store exon column indices
		
		for i in "${!header_cols[@]}"; do
		    case "${header_cols[$i]}" in
			"contig") contig_col=$i ;;
			C-EXON_[1-9]_start) exon_cols["${header_cols[$i]}"]=$i ;;  # Match any C-EXON_X_start where X is 1-9
			C-EXON_[1-9]_end) exon_cols["${header_cols[$i]}"]=$i ;;  # Match any C-EXON_X_end where X is 1-9
		    esac
		done
		
		# To access a specific exon start or end column index later in the script, you can use:
		# ${exon_cols[C-EXON_4_start]} or ${exon_cols[C-EXON_4_end]}, etc.
		
		# Specify the output file name
		output_file="exon_regions.tsv"
		
		# Write the header to the output file
	#	echo -e "Contig\tExon_Start\tExon_End" 
		> "$output_file"
		
		# Iterate over each line of the input file, skipping the header
		tail -n +2 "$modified_import_out" | while IFS=, read -ra line; do
		    contig="${line[$contig_col]}"
		    
		    # Generate a unique filename for the BED file based on the contig name
		    contig_filename=$(echo "$contig" | tr '/' '_')
		    bed_file="${base_outd}/${gene_type}/${contig_filename}.bed"
		    mkdir -p "$(dirname "$bed_file")"
		    
		    # Clear the BED file for new entries
		    > "$bed_file"
		    
		    # Iterate over each exon start and end index pair from the exon_cols associative array
		    for key in "${!exon_cols[@]}"; do
			if [[ $key =~ C-EXON_[1-9]_start ]]; then
			    exon_start=$(echo "${line[${exon_cols[$key]}]}" | awk '{printf "%.0f", $1}')
			    exon_end_key="${key/_start/_end}"  # Replace 'start' with 'end' to get the corresponding end key
			    exon_end=$(echo "${line[${exon_cols[$exon_end_key]}]}" | awk '{printf "%.0f", $1}')
			    
			    # Write the contig, exon start, and exon end to the BED file
			    echo -e "${contig}\t${exon_start}\t${exon_end}" >> "$bed_file"
			fi
		    done
		    

#                    contig_filename=$(echo "$contig" | tr '/' '_')
                    tmp_bam="${base_outd}/${gene_type}/${contig_filename}_${start}_${end}.bam"
                    mkdir -p "$(dirname "$tmp_bam")"
                    samtools view -F 0x100 -F 0x800 -b "$bam_file" -o "$tmp_bam" -U "/dev/null" "${contig}:${start}-${end}"
                    samtools index "$tmp_bam"
                    
		    total_positions=$(awk '{sum += $3 - $2} END {print sum}' "$bed_file")
                    samtools mpileup -f "$ref" -l "$bed_file" "$tmp_bam" | \
			awk -v total_positions="$total_positions" \
                    'BEGIN {
                        total_reads=0; mismatched_positions=0; matched_positions=0;
                        mismatch_list=""; match_list="";
                        mismatched_positions_coverage_less_than_10=0;
                        mismatched_positions_coverage_10_or_greater=0;
                        matched_positions_coverage_less_than_10=0;
                        matched_positions_coverage_10_or_greater=0;
                    }
                    {
                        total_reads += length($5);
                        mismatches = length(gensub(/[.,]/, "", "g", $5));
                        matches = length(gensub(/[^.,]/, "", "g", $5));
                        mismatch_list = (mismatch_list == "" ? mismatches : mismatch_list ":" mismatches);
                        match_list = (match_list == "" ? matches : match_list ":" matches);

                        coverage = length($5);
                        mismatch_rate = mismatches / coverage;
                        match_rate = matches / coverage;

                        if (mismatch_rate > 0.2) {
                            mismatched_positions++;
                            if (coverage < 10) {
                                mismatched_positions_coverage_less_than_10++;
                            } else {
                                mismatched_positions_coverage_10_or_greater++;
                            }
                        }

                        if (match_rate > 0.8) {
                            matched_positions++;
                            if (coverage < 10) {
                                matched_positions_coverage_less_than_10++;
                            } else {
                                matched_positions_coverage_10_or_greater++;
                            }
                        }
                    }
                    END {
                        percent_accuracy = (matched_positions / total_positions) * 100;
                        print total_positions, total_reads, mismatched_positions, matched_positions, mismatch_list, match_list, mismatched_positions_coverage_less_than_10, mismatched_positions_coverage_10_or_greater, matched_positions_coverage_less_than_10, matched_positions_coverage_10_or_greater, percent_accuracy;
                    }' OFS=',' >> "$tmp_file"

                    rm "$tmp_bam" "${tmp_bam}.bai"
                done

                paste -d ',' "$modified_import_out" "$tmp_file" > "${modified_import_out%.csv}_with_read_support.csv"
                rm -f "$tmp_file" "$modified_import_out"
            fi
        done
    done < extended_samples_paths.fofn
}


#run_make_ref_masked
#run_map_ccs_to_pers
#run_append_pos
get_read_support_vdj3
# do not use this get_read_support_ighc
