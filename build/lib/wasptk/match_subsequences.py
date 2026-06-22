"""Utilities for per-read sequence matching."""
from typing import Tuple, Optional
import pandas as pd
import pysam

_complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

def reverse_complement(seq: str) -> str:
    return "".join(_complement.get(b, "N") for b in seq[::-1])


def extract_sequence(
    row: pd.Series,
    gene: str,
    vseq_col: Optional[str] = None,
    dseq_col: Optional[str] = None,
    jseq_col: Optional[str] = None,
    cseq_col: Optional[str] = None,
) -> str:
    """Return a sequence for ``gene`` from ``row`` respecting strand.

    The fourth character of ``gene`` is used to pick a column from the provided
    names. Missing columns or values result in an empty string.
    """

    type_map = {"V": vseq_col, "D": dseq_col, "J": jseq_col, "C": cseq_col}
    gene_type = gene[3].upper() if len(gene) >= 4 else ""
    seq_col = type_map.get(gene_type)
    if not seq_col:
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
