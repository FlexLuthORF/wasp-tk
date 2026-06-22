"""Compute read support metrics from a mapped BAM.

``compute_read_support`` reads an annotation table and calculates coverage and
per-read sequence matching statistics. Sequence columns can be specified
separately for V, D, J and C genes. These columns default to ``gene_seq``.
"""

from __future__ import annotations

from typing import List, Dict, Tuple, Optional
import subprocess

import pandas as pd
import pysam

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement.get(base, 'N') for base in seq[::-1]])

def get_sequence_from_row(row, gene_key):
    """Extracts sequence string from CSV row."""
    reverse_comp = False
    if 'sense' in row and row['sense'] == '-':
        reverse_comp = True

    sequence_column = None
    if any(s in gene_key for s in ['IGKV', 'IGLV', 'IGHV', 'TRAV', 'TRBV', 'TRDV', 'TRGV']):
        sequence_column = 'V-REGION'
    elif any(s in gene_key for s in ['IGKJ', 'IGLJ', 'IGHJ', 'TRAJ', 'TRBJ', 'TRDJ', 'TRGJ']):
        sequence_column = 'J-REGION'
    elif any(s in gene_key for s in ['IGHD', 'TRBD', 'TRDD']):
        sequence_column = 'D-REGION'
    elif any(s in gene_key for s in ['IGKC', 'IGLC', 'TRAC', 'TRBC', 'TRDC', 'TRGC']):
        sequence_column = 'C-REGION'
    
    if not sequence_column or sequence_column not in row:
        return ""

    sequence = row[sequence_column]
    if reverse_comp and sequence:
        sequence = reverse_complement(sequence)
    return sequence

def count_matching_reads_internal(bamfile, chrom, start, end, sequence):
    """Counts full span and perfect matches using pysam."""
    full_span_count = 0
    perfect_match_count = 0
    try:
        for read in bamfile.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_start <= start and read.reference_end >= end:
                full_span_count += 1
                read_seq = read.query_sequence
                if sequence and (sequence in read_seq or reverse_complement(sequence) in read_seq):
                    perfect_match_count += 1
    except Exception as e:
        print(f"[WARNING] Pysam error in count_matching_reads: {e}")
        return 0, 0
    return full_span_count, perfect_match_count

def get_ighc_sequences_from_row(row):
    """Extracts IGHC exon sequences."""
    sequences = []
    exon_start_end_list = []
    reverse_comp = 'sense' in row and row['sense'] == '-'
    
    for i in range(1, 10):
        seq_col = f'C-EXON_{i}'
        start_col = f'C-EXON_{i}_start'
        end_col = f'C-EXON_{i}_end'
        sequence = row.get(seq_col, '')
        start_val = row.get(start_col, '')
        end_val = row.get(end_col, '')

        if start_val and end_val:
            try:
                s_int = int(float(start_val))
                e_int = int(float(end_val))
                exon_start_end_list.append((s_int, e_int))
            except ValueError: pass
        
        if reverse_comp and sequence:
            sequence = reverse_complement(sequence)
        if sequence:
            sequences.append(sequence)
            
    return sequences, exon_start_end_list

def count_ighc_matches_internal(bamfile, contig, start, end, sequences, exon_start_end_list):
    """IGHC matching logic."""
    full_span_count = 0
    full_span_all_match_count = 0
    
    max_len = max(len(sequences), len(exon_start_end_list))
    if max_len == 0: return 0, 0, [], []
    
    perfect_matches = [0] * len(sequences)
    perfect_spans = [0] * len(exon_start_end_list)

    try:
        for read in bamfile.fetch(contig, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
            
            if read.reference_start <= start and read.reference_end >= end:
                full_span_count += 1
                read_seq = read.query_sequence
                all_matches_flag = True
                for seq in sequences:
                    if not seq: continue
                    if (seq not in read_seq) and (reverse_complement(seq) not in read_seq):
                        all_matches_flag = False
                        break
                if all_matches_flag and sequences:
                    full_span_all_match_count += 1
            
            read_seq = read.query_sequence
            for idx, (ex_s, ex_e) in enumerate(exon_start_end_list):
                if read.reference_start <= ex_s and read.reference_end >= ex_e:
                    perfect_spans[idx] += 1
            for idx, seq in enumerate(sequences):
                if seq and (seq in read_seq or reverse_complement(seq) in read_seq):
                    perfect_matches[idx] += 1
    except Exception as e:
        print(f"[WARNING] Pysam error in ighc count: {e}")

    return full_span_count, full_span_all_match_count, perfect_matches, perfect_spans

_CONSTANT_GENES = ["IGKC", "IGLC", "TRAC", "TRBC", "TRDC", "TRGC", "IGHC"]

def validate_vdj_stats(sample_name, gene_name, stats_tuple):
    """
    Sanity Check: Validates that the summary statistics (columns) match the 
    raw colon-separated position data strings.
    """
    expected_total_positions = stats_tuple[0]
    reported_mismatch_count  = stats_tuple[2]
    reported_match_count     = stats_tuple[3]
    reported_10x_count       = stats_tuple[7]
    
    raw_mismatch_string = stats_tuple[4]
    raw_match_string    = stats_tuple[5]
    
    mismatch_list = raw_mismatch_string.split(':') if raw_mismatch_string else []
    match_list    = raw_match_string.split(':')    if raw_match_string    else []

    if len(mismatch_list) != expected_total_positions or len(match_list) != expected_total_positions:
        print(f"[SANITY FAIL] {sample_name} {gene_name}: Data gap detected!")
        print(f"  Expected {expected_total_positions} positions based on coordinates.")
        print(f"  Found {len(match_list)} match entries and {len(mismatch_list)} mismatch entries.")
        return False

    calculated_match_count = 0
    calculated_mismatch_count = 0
    calculated_10x_count = 0

    for match_str_val, mismatch_str_val in zip(match_list, mismatch_list):
        try:
            match_read_count = int(float(match_str_val))
            mismatch_read_count = int(float(mismatch_str_val))
            
            total_coverage_at_pos = match_read_count + mismatch_read_count
            
            if total_coverage_at_pos > 0:
                match_ratio = match_read_count / total_coverage_at_pos
                if match_ratio > 0.8:
                    calculated_match_count += 1
                
                mismatch_ratio = mismatch_read_count / total_coverage_at_pos
                if mismatch_ratio > 0.2:
                    calculated_mismatch_count += 1
            
            if total_coverage_at_pos >= 10:
                calculated_10x_count += 1
                
        except ValueError:
            continue

    errors = []
    
    if reported_match_count != calculated_match_count:
        errors.append(f"Matched_Positions mismatch: Header says {reported_match_count}, Raw Data says {calculated_match_count}")
        
    if reported_mismatch_count != calculated_mismatch_count:
        errors.append(f"Mismatched_Positions mismatch: Header says {reported_mismatch_count}, Raw Data says {calculated_mismatch_count}")
        
    if reported_10x_count != calculated_10x_count:
        errors.append(f"10x_Coverage mismatch: Header says {reported_10x_count}, Raw Data says {calculated_10x_count}")

    if errors:
        print(f"[SANITY FAIL] {sample_name} {gene_name} Logic Error:")
        for error in errors:
            print(f"  - {error}")
        return False

    return True

def _pileup_region_lines(
    bam_path: str,
    contig: str,
    start: int,
    end: int,
    reference: str,
) -> List[str]:
    cmd = ["samtools", "mpileup", "-f", reference]
    cmd += ["-r", f"{contig}:{start}-{end}", bam_path]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[WARNING] mpileup failed: {res.stderr}")
    return [line for line in res.stdout.split('\n') if line.strip()]

def parse_mpileup_and_calculate(lines: List[str], total_positions: int):
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    positions_with_10x = 0
    mismatch_list = []
    match_list = []

    for line in lines:
        if not line.strip(): continue
        parts = line.strip().split()
        if len(parts) < 5: continue
        coverage = int(parts[3])
        read_bases = parts[4]
        total_reads += coverage

        match_count = sum(1 for c in read_bases if c in ['.', ','])
        mismatch_count = coverage - match_count
        mismatch_list.append(str(mismatch_count))
        match_list.append(str(match_count))

        if coverage >= 10: positions_with_10x += 1
        
        mismatch_rate = mismatch_count / coverage if coverage > 0 else 0
        match_rate = match_count / coverage if coverage > 0 else 0
        
        if mismatch_rate > 0.2: mismatched_positions += 1
        if match_rate > 0.8: matched_positions += 1

    avg_reads_per_position = total_reads / total_positions if total_positions > 0 else 0
    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0

    return (total_positions, avg_reads_per_position, mismatched_positions, matched_positions,
            ":".join(mismatch_list), ":".join(match_list), percent_accuracy, positions_with_10x)

def parse_mpileup_ighc(lines: List[str], total_positions: int):
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    mismatch_list = []
    match_list = []
    mism_lt10 = mism_ge10 = mat_lt10 = mat_ge10 = 0

    for line in lines:
        if not line.strip(): continue
        parts = line.strip().split()
        if len(parts) < 5: continue
        coverage = int(parts[3])
        read_bases = parts[4]
        total_reads += coverage
        
        matches = sum(1 for c in read_bases if c in ['.', ','])
        mismatches = coverage - matches
        mismatch_list.append(str(mismatches))
        match_list.append(str(matches)) 

        mismatch_rate = mismatches / coverage if coverage > 0 else 0
        match_rate = matches / coverage if coverage > 0 else 0

        if mismatch_rate > 0.2:
            mismatched_positions += 1
            if coverage < 10: mism_lt10 += 1
            else: mism_ge10 += 1

        if match_rate >= 0.8:
            matched_positions += 1
            if coverage < 10: mat_lt10 += 1
            else: mat_ge10 += 1

    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0
    avg_cov = (total_reads/total_positions) if total_positions > 0 else 0
    pos_10x = (mism_ge10 + mat_ge10)

    return (total_positions, total_reads, avg_cov, mismatched_positions, matched_positions,
            ":".join(mismatch_list), ":".join(match_list), mism_lt10, mism_ge10, mat_lt10, mat_ge10,
            pos_10x, percent_accuracy)

def compute_read_support(
    allele_table: str,
    bam_path: str,
    output: str,
    reference: str,
    contig_col: str = "contig",
    start_col: str = "start",
    end_col: str = "end",
    gene_col: str = "gene",
    vseq_col: str = "gene_seq",
    dseq_col: str = "gene_seq",
    jseq_col: str = "gene_seq",
    cseq_col: str = "gene_seq",
) -> None:
    df = pd.read_csv(allele_table)
    bam = pysam.AlignmentFile(bam_path, "rb")

    metrics: List[Dict[str, int | float | str]] = []

    for _, row in df.iterrows():
        contig = str(row[contig_col])
        start = int(row[start_col])
        end = int(row[end_col])
        gene = str(row.get(gene_col, ""))
        sample_name = str(row.get('sample_name', 'Unknown'))

        if any(gene.startswith(c) for c in _CONSTANT_GENES):
            seqs, coords = get_ighc_sequences_from_row(row)
            all_lines: List[str] = []
            total_positions = 0
            if coords:
                for s, e in coords:
                    lines = _pileup_region_lines(bam_path, contig, s, e, reference)
                    all_lines.extend(lines)
                    total_positions += (e - s + 1)
            else:
                all_lines = _pileup_region_lines(bam_path, contig, start, end, reference)
                total_positions = end - start + 1
            (
                tpos,
                total_reads,
                avg_cov,
                mismatched_positions,
                matched_positions,
                mismatch_str,
                match_str,
                mismatched_positions_cov_lt10,
                mismatched_positions_cov_ge10,
                matched_positions_cov_lt10,
                matched_positions_cov_ge10,
                pos_10x,
                pct_acc,
            ) = parse_mpileup_ighc(all_lines, total_positions)
            full_span, all_match, match_counts, span_counts = count_ighc_matches_internal(
                bam, contig, start, end, seqs, coords
            )
            row_metrics: Dict[str, int | float | str] = {
                "Total_Positions": tpos,
                "Total_Reads_by_Positions": total_reads,
                "Average_Coverage": avg_cov,
                "Mismatched_Positions": mismatched_positions,
                "Matched_Positions": matched_positions,
                "Position_Mismatches": mismatch_str,
                "Position_Matches": match_str,
                "Mismatched_Positions_Coverage_Less_Than_10": mismatched_positions_cov_lt10,
                "Mismatched_Positions_Coverage_10_Or_Greater": mismatched_positions_cov_ge10,
                "Matched_Positions_Coverage_Less_Than_10": matched_positions_cov_lt10,
                "Matched_Positions_Coverage_10_Or_Greater": matched_positions_cov_ge10,
                "Positions_With_At_Least_10x_Coverage": pos_10x,
                "Percent_Accuracy": pct_acc,
                "Fully_Spanning_Reads": full_span,
                "Fully_Spanning_Reads_100%_Match": all_match,
            }
            for i, count in enumerate(match_counts, start=1):
                row_metrics[f"Allele_reads_100_Match_e{i}"] = count
            for i, count in enumerate(span_counts, start=1):
                row_metrics[f"Allele_reads_fully_spanning_e{i}"] = count
        else:
            total_positions = end - start + 1
            lines = _pileup_region_lines(bam_path, contig, start, end, reference)
            stats = parse_mpileup_and_calculate(lines, total_positions)
            
            validate_vdj_stats(sample_name, gene, stats)
            
            (
                tpos,
                avg_cov,
                mismatched_positions,
                matched_positions,
                mismatch_str,
                match_str,
                pct_acc,
                pos_10x,
            ) = stats
            
            sequence = get_sequence_from_row(row, gene)
            
            if not sequence or len(sequence) < 1:
                print(f"[WARNING] {sample_name} {gene}: Sequence missing or short (<1bp).")
                
            if sequence:
                full_span, perfect = count_matching_reads_internal(bam, contig, start, end, sequence)
            else:
                full_span, perfect = 0, 0
                
            row_metrics = {
                "Total_Positions": tpos,
                "Average_Coverage": avg_cov,
                "Mismatched_Positions": mismatched_positions,
                "Matched_Positions": matched_positions,
                "Position_Mismatches": mismatch_str,
                "Position_Matches": match_str,
                "Percent_Accuracy": pct_acc,
                "Positions_With_At_Least_10x_Coverage": pos_10x,
                "Fully_Spanning_Reads": full_span,
                "Fully_Spanning_Reads_100%_Match": perfect,
            }

        metrics.append(row_metrics)

    all_cols: List[str] = []
    for m in metrics:
        for k in m.keys():
            if k not in all_cols:
                all_cols.append(k)

    for col in all_cols:
        df[col] = [m.get(col, 0) for m in metrics]

    df.to_csv(output, index=False)
