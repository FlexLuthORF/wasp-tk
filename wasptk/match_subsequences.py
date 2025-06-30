"""Utilities for per-read sequence matching."""
from typing import Tuple, Optional
import pandas as pd
import pysam

_complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

def reverse_complement(seq: str) -> str:
    return "".join(_complement.get(b, "N") for b in seq[::-1])


def extract_sequence(row: pd.Series, seq_col: Optional[str]) -> str:
    """Return the sequence from ``seq_col`` respecting strand information.

    ``seq_col`` may be ``None`` to disable extraction.
    """

    if seq_col is None:
        return ""

    if pd.isna(row.get("sense")):
        rev = False
    else:
        rev = str(row.get("sense")) == "-"

    if seq_col not in row or pd.isna(row[seq_col]):
        return ""

    seq = str(row[seq_col])
    if rev:
        seq = reverse_complement(seq)
    return seq


def count_matching_reads(bam: pysam.AlignmentFile, chrom: str, start: int, end: int, sequence: str) -> Tuple[int, int]:
    full_span = 0
    perfect = 0
    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start <= start and read.reference_end >= end:
            full_span += 1
            seq = read.query_sequence
            if sequence and (sequence in seq or reverse_complement(sequence) in seq):
                perfect += 1
    return full_span, perfect
