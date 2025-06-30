import argparse
import json
from .read_support import compute_read_support
from .ancestry import run_aims


def _cmd_readsupport(args: argparse.Namespace) -> None:
    compute_read_support(
        args.allele_table,
        args.bam,
        args.output,
        reference=args.reference,
        contig_col=args.contig_col,
        start_col=args.start_col,
        end_col=args.end_col,
        gene_col=args.gene_col,
    )


def _cmd_plotcov(args: argparse.Namespace) -> None:
    from .plotcov import plotcov
    plotcov(args.depth, args.loci, args.out)


def _cmd_aims(args: argparse.Namespace) -> None:
    res = run_aims(args.vcf, sample_id=args.sample, bam=args.bam)
    if args.output:
        with open(args.output, "w") as fh:
            json.dump(res, fh, indent=2)
    else:
        print(json.dumps(res, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(prog="wasptk", description="WASP toolkit")
    sub = parser.add_subparsers(dest="command", required=True)

    p_read = sub.add_parser("readsupport", help="Compute read support from a BAM")
    p_read.add_argument("allele_table", help="Path to allele_annotation.csv")
    p_read.add_argument("bam", help="Mapped reads BAM")
    p_read.add_argument("output", help="Output CSV")
    p_read.add_argument("--reference", help="Reference FASTA for samtools mpileup")
    p_read.add_argument("--contig-col", default="contig", help="Column for contig name")
    p_read.add_argument("--start-col", default="start", help="Column for start position")
    p_read.add_argument("--end-col", default="end", help="Column for end position")
    p_read.add_argument("--gene-col", default="gene", help="Column for gene name")
    p_read.set_defaults(func=_cmd_readsupport)

    p_plot = sub.add_parser(
        "plotcov", help="Plot coverage graphs from mosdepth output"
    )
    p_plot.add_argument("--depth", required=True, help="mosdepth per-base BED.gz")
    p_plot.add_argument(
        "--loci", required=True, help="4-column BED: chrom, start, end, name"
    )
    p_plot.add_argument(
        "--out",
        required=True,
        help="output images prefix (e.g. PNG)",
    )
    p_plot.set_defaults(func=_cmd_plotcov)

    p_aims = sub.add_parser("aims", help="Infer ancestry using AIMs")
    p_aims.add_argument(
        "-v",
        "--vcf",
        required=True,
        help="Input VCF with AIM variants",
    )
    p_aims.add_argument(
        "-b",
        "--bam",
        help="BAM file for coverage calculation",
    )
    p_aims.add_argument(
        "-s",
        "--sample",
        help="Sample identifier to include in the output",
    )
    p_aims.add_argument(
        "-o",
        "--output",
        help="Write JSON result to file instead of stdout",
    )
    p_aims.set_defaults(func=_cmd_aims)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
