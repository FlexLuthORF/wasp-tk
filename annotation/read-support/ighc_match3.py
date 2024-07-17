#!/usr/bin/env python3

import csv
import sys
import pysam

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])

def get_sequences_and_regions_from_csv(import_csv, gene_key, contig):
    sequences = []
    regions = []
    with open(import_csv, mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            if row['gene'] == gene_key and row['contig'] == contig:
                start = int(float(row['allele_sequence_start']))
                end = int(float(row['allele_sequence_end']))
                reverse_comp = 'sense' in row and row['sense'] == '-'
                # Populate sequences from C-EXON columns
                for i in range(1, 10):
                    sequence_column = f'C-EXON_{i}'
                    sequence = row.get(sequence_column, '')
                    if reverse_comp and sequence:
                        sequence = reverse_complement(sequence)
                    if sequence:  # Only add non-empty sequences
                        sequences.append(sequence)
                regions.append((contig, start, end))  # Make sure to append a tuple of (contig, start, end)
                break
    return sequences, regions

def count_matching_reads(bamfile, sequences, regions):
    full_span_counts = []
    perfect_match_counts = []
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for region in regions:
        start, end = region[1], region[2]  # Assuming region tuple structure is (contig, start, end)
        contig = str(region[0])  # Ensure this is a string as expected
        full_span_count = 0
        perfect_matches = [0] * len(sequences)
        for read in samfile.fetch(str(contig), start, end):  # Pass contig correctly
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_start <= start and read.reference_end >= end:
                full_span_count += 1
                read_seq = read.query_sequence
                for idx, sequence in enumerate(sequences):
                    if sequence and (sequence in read_seq or reverse_complement(sequence) in read_seq):
                        perfect_matches[idx] += 1
        full_span_counts.append(full_span_count)
        perfect_match_counts.append(perfect_matches)
    return full_span_counts, perfect_match_counts

if __name__ == '__main__':
    bamfile = sys.argv[1]
    contig = sys.argv[2]
    gene_key = sys.argv[3]
    import_csv = sys.argv[4]

    sequences, regions = get_sequences_and_regions_from_csv(import_csv, gene_key, contig)
    if sequences and regions:
        full_span_counts, perfect_match_counts = count_matching_reads(bamfile, sequences, regions)
        for i, counts in enumerate(perfect_match_counts):
            print(f"{full_span_counts[i]},{','.join(map(str, counts))}")
    else:
        print(f"No sequences or regions found for gene {gene_key} on contig {contig}", file=sys.stderr)
