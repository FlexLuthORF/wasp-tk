# WASP Toolkit

This repository provides scripts and utilities for immune receptor analysis. A Python
package exposes a `wasptk-read-support` command that computes read support metrics
from a mapped BAM and an annotation table.

## Installation

```bash
pip install -e .
```

## Usage

```bash
wasptk-read-support <allele_annotation.csv> <mapped.bam> <reference.fasta> <output.csv>
```

Optional flags allow overriding the column names in the annotation table using
`--contig-col`, `--start-col`, `--end-col`, and `--gene-col`.
