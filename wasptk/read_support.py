"""Compute read support metrics from a mapped BAM."""

from __future__ import annotations

from typing import List, Dict, Tuple
import os
import subprocess
import tempfile

import pandas as pd
import pysam

from .match_subsequences import extract_sequence, count_matching_reads as count_match_sub
from .ighc_match import extract_exon_sequences, count_matching_reads as count_ighc

_CONSTANT_GENES = ["IGKC", "IGLC", "TRAC", "TRBC", "TRDC", "TRGC", "IGHC"]


def run_mpileup(bam: str, reference: str, region: str | None = None, bed: str | None = None) -> List[str]:
    cmd = ["samtools", "mpileup", "-f", reference]
    if bed:
        cmd += ["-l", bed]
    elif region:
        cmd += ["-r", region]
    cmd.append(bam)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr)
    return result.stdout.strip().splitlines()


def parse_mpileup_and_calculate(lines: List[str], total_positions: int) -> Tuple[int, float, int, int, str, str, float, int]:
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    positions_with_10x = 0
    mismatch_list: List[str] = []
    match_list: List[str] = []

    for line in lines:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        read_bases = parts[4]
        coverage = len(read_bases)
        total_reads += coverage

        mismatch_count = sum(1 for c in read_bases if c not in [".", ","])
        match_count = coverage - mismatch_count

        mismatch_list.append(str(mismatch_count))
        match_list.append(str(match_count))

        if coverage >= 10:
            positions_with_10x += 1

        mismatch_rate = mismatch_count / coverage if coverage > 0 else 0
        match_rate = match_count / coverage if coverage > 0 else 0

        if mismatch_rate > 0.2:
            mismatched_positions += 1
        if match_rate > 0.8:
            matched_positions += 1

    avg_reads_per_position = total_reads / total_positions if total_positions > 0 else 0
    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0

    mismatch_str = ":".join(mismatch_list)
    match_str = ":".join(match_list)

    return (
        total_positions,
        avg_reads_per_position,
        mismatched_positions,
        matched_positions,
        mismatch_str,
        match_str,
        percent_accuracy,
        positions_with_10x,
    )


def parse_mpileup_ighc(lines: List[str], total_positions: int) -> Tuple[int, int, float, int, int, str, str, int, int, int, int, int, float]:
    total_reads = 0
    mismatched_positions = 0
    matched_positions = 0
    mismatch_list: List[str] = []
    match_list: List[str] = []
    mismatched_positions_cov_lt10 = 0
    mismatched_positions_cov_ge10 = 0
    matched_positions_cov_lt10 = 0
    matched_positions_cov_ge10 = 0

    for line in lines:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) < 5:
            continue
        coverage = int(parts[3])
        read_bases = parts[4]
        total_reads += coverage

        matches = sum(1 for c in read_bases if c in [".", ","])
        mismatches = coverage - matches

        mismatch_list.append(str(mismatches))
        match_list.append(str(matches))

        mismatch_rate = mismatches / coverage if coverage > 0 else 0
        match_rate = matches / coverage if coverage > 0 else 0

        if mismatch_rate > 0.2:
            mismatched_positions += 1
            if coverage < 10:
                mismatched_positions_cov_lt10 += 1
            else:
                mismatched_positions_cov_ge10 += 1

        if match_rate >= 0.8:
            matched_positions += 1
            if coverage < 10:
                matched_positions_cov_lt10 += 1
            else:
                matched_positions_cov_ge10 += 1

    percent_accuracy = (matched_positions / total_positions) * 100 if total_positions > 0 else 0
    avg_cov = (total_reads / total_positions) if total_positions > 0 else 0
    pos_10x = mismatched_positions_cov_ge10 + matched_positions_cov_ge10

    return (
        total_positions,
        total_reads,
        avg_cov,
        mismatched_positions,
        matched_positions,
        ":".join(mismatch_list),
        ":".join(match_list),
        mismatched_positions_cov_lt10,
        mismatched_positions_cov_ge10,
        matched_positions_cov_lt10,
        matched_positions_cov_ge10,
        pos_10x,
        percent_accuracy,
    )


def compute_read_support(
    allele_table: str,
    bam_path: str,
    reference: str,
    output: str,
    contig_col: str = "contig",
    start_col: str = "start",
    end_col: str = "end",
    gene_col: str = "gene",
) -> None:
    df = pd.read_csv(allele_table)
    bam = pysam.AlignmentFile(bam_path, "rb")

    metrics: List[Dict[str, int | float | str]] = []

    for _, row in df.iterrows():
        contig = str(row[contig_col])
        start = int(row[start_col])
        end = int(row[end_col])
        gene = str(row.get(gene_col, ""))
        region = f"{contig}:{start}-{end}"

        if any(gene.startswith(c) for c in _CONSTANT_GENES):
            seqs, coords = extract_exon_sequences(row)
            if coords:
                with tempfile.NamedTemporaryFile(mode="w", delete=False) as bed:
                    for s, e in coords:
                        bed.write(f"{contig}\t{s-1}\t{e}\n")
                    bed_path = bed.name
                lines = run_mpileup(bam_path, reference, bed=bed_path)
                os.unlink(bed_path)
                total_pos = sum(e - (s - 1) for s, e in coords)
            else:
                lines = run_mpileup(bam_path, reference, region)
                total_pos = end - start + 1
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
            ) = parse_mpileup_ighc(lines, total_pos)
            full_span, all_match, match_counts, span_counts = count_ighc(
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
            lines = run_mpileup(bam_path, reference, region)
            total_pos = end - start + 1
            (
                tpos,
                avg_cov,
                mismatched_positions,
                matched_positions,
                mismatch_str,
                match_str,
                pct_acc,
                pos_10x,
            ) = parse_mpileup_and_calculate(lines, total_pos)
            seq = extract_sequence(row, gene)
            full_span, perfect = count_match_sub(bam, contig, start, end, seq)
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
