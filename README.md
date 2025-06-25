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
