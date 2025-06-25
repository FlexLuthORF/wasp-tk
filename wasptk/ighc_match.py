"""Subset of functionality from annotation/read-support/ighc_match3.py."""

from typing import List, Tuple

import pandas as pd
import pysam

from .match_subsequences import reverse_complement


def extract_exon_sequences(row: pd.Series) -> Tuple[List[str], List[Tuple[int, int]]]:
    """Return lists of exon sequences and coordinates from a row."""
    seqs = []
    coords = []
    for i in range(1, 10):
        col = f"C-EXON_{i}"
        s_col = f"C-EXON_{i}_start"
        e_col = f"C-EXON_{i}_end"
        if col in row and pd.notna(row[col]) and row[col] != "":
            seq = str(row[col])
            if row.get("sense") == "-":
                seq = reverse_complement(seq)
            seqs.append(seq)
            if s_col in row and e_col in row and pd.notna(row[s_col]) and pd.notna(row[e_col]):
                coords.append((int(row[s_col]), int(row[e_col])))
    return seqs, coords


def count_matching_reads(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    sequences: List[str],
    exon_coords: List[Tuple[int, int]],
) -> Tuple[int, int]:
    """Return counts of reads spanning region and matching all exon sequences."""
    full_span = 0
    full_span_all_match = 0
    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start <= start and read.reference_end >= end:
            full_span += 1
            seq = read.query_sequence
            if all(seq_part in seq or reverse_complement(seq_part) in seq for seq_part in sequences):
                full_span_all_match += 1
    return full_span, full_span_all_match
