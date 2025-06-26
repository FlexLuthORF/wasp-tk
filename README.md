# WASP Toolkit

This repository provides a lightweight Python CLI package for computing read support metrics.

## Installation

```bash
pip install -e .
```

## Usage

Run the `wasptk-read-support` command on an allele annotation table and a BAM file of mapped reads:

```bash
wasptk-read-support <allele_annotation.csv> <mapped.bam> <output.csv>
```

Optional flags allow overriding the column names used from the annotation table:
`--contig-col`, `--start-col`, `--end-col`, and `--gene-col`.

The output table includes two columns for most genes: `full_span_reads` and `perfect_match_reads`. For constant-region genes (e.g. IGKC/IGLC/TRAC/TRBC/TRDC/TRGC) additional columns report per-exon support, such as `all_exons_match_reads`, `exon1_match_reads`, and `exon1_span_reads`.
