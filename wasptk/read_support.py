"""Compute read support metrics from a mapped BAM."""

from typing import List

import pandas as pd
import pysam

from .match_subsequences import extract_sequence, count_matching_reads as count_match_sub
from .ighc_match import extract_exon_sequences, count_matching_reads as count_ighc


_CONSTANT_GENES = ["IGKC", "IGLC", "TRAC", "TRBC", "TRDC", "TRGC"]


def compute_read_support(
    allele_table: str,
    bam_path: str,
    output: str,
    contig_col: str = "contig",
    start_col: str = "start",
    end_col: str = "end",
    gene_col: str = "gene",
) -> None:
    """Add read support columns to ``allele_table`` and write ``output``."""
    df = pd.read_csv(allele_table)
    bam = pysam.AlignmentFile(bam_path, "rb")

    full_span_list: List[int] = []
    perfect_list: List[int] = []

    for _, row in df.iterrows():
        contig = str(row[contig_col])
        start = int(row[start_col])
        end = int(row[end_col])
        gene = str(row.get(gene_col, ""))

        if any(gene.startswith(c) for c in _CONSTANT_GENES):
            sequences, _coords = extract_exon_sequences(row)
            full_span, perfect = count_ighc(bam, contig, start, end, sequences, _coords)
        else:
            seq = extract_sequence(row, gene)
            full_span, perfect = count_match_sub(bam, contig, start, end, seq)

        full_span_list.append(full_span)
        perfect_list.append(perfect)

    df["full_span_reads"] = full_span_list
    df["perfect_match_reads"] = perfect_list
    df.to_csv(output, index=False)
