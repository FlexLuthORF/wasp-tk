"""Compute read support metrics from a mapped BAM."""

from typing import List, Dict

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

    metrics: List[Dict[str, int]] = []

    for _, row in df.iterrows():
        contig = str(row[contig_col])
        start = int(row[start_col])
        end = int(row[end_col])
        gene = str(row.get(gene_col, ""))

        if any(gene.startswith(c) for c in _CONSTANT_GENES):
            sequences, coords = extract_exon_sequences(row)
            full_span, all_match, match_counts, span_counts = count_ighc(
                bam, contig, start, end, sequences, coords
            )
            row_metrics: Dict[str, int] = {
                "full_span_reads": full_span,
                "all_exons_match_reads": all_match,
            }
            for i, count in enumerate(match_counts, start=1):
                row_metrics[f"exon{i}_match_reads"] = count
            for i, count in enumerate(span_counts, start=1):
                row_metrics[f"exon{i}_span_reads"] = count
        else:
            seq = extract_sequence(row, gene)
            full_span, perfect = count_match_sub(bam, contig, start, end, seq)
            row_metrics = {
                "full_span_reads": full_span,
                "perfect_match_reads": perfect,
            }

        metrics.append(row_metrics)

    # collect all metric keys to ensure consistent columns
    all_cols = set()
    for m in metrics:
        all_cols.update(m.keys())
    for col in sorted(all_cols):
        df[col] = [m.get(col, 0) for m in metrics]
    df.to_csv(output, index=False)
