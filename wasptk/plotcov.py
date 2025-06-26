"""Coverage plotting utilities for mosdepth output."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import ScalarFormatter


# ----- shared helpers -----

def load_loci(path: Path) -> pd.DataFrame:
    """Read a 4-column BED (chrom, start, end, name) sorted by chrom/start."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
        dtype={"chrom": str},
    )
    df.sort_values(["chrom", "start"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


# ----- depth facet plot -----

def collect_depth(depth_gz: Path, loci_df: pd.DataFrame) -> List[List[Tuple[int, int]]]:
    """Return per-locus depth coordinates for facet plotting."""
    facet_data: List[List[Tuple[int, int]]] = [[] for _ in range(len(loci_df))]
    by_chrom: dict[str, List[Tuple[int, int, int]]] = {}
    for idx, row in loci_df.iterrows():
        by_chrom.setdefault(row.chrom, []).append((idx, row.start, row.end))

    with gzip.open(depth_gz, "rt") as fh:
        for ln in fh:
            chrom, s_str, e_str, d_str = ln.strip().split("\t")
            s = int(s_str)
            e = int(e_str)
            d = int(d_str)
            for idx, L, R in by_chrom.get(chrom, []):
                if e > L and s < R:
                    facet_data[idx].append((s, d))
                    facet_data[idx].append((e, d))

    for lst in facet_data:
        lst[:] = sorted(set(lst), key=lambda t: t[0])
    return facet_data


def plot_depth_facets(
    depth_gz: Path,
    loci_df: pd.DataFrame,
    depth_data: List[List[Tuple[int, int]]],
    out_png: Path,
) -> None:
    n = len(loci_df)
    fig, axes = plt.subplots(
        nrows=n,
        ncols=1,
        figsize=(12, 1.5 * n),
        sharey=True,
        constrained_layout=True,
    )
    if n == 1:
        axes = [axes]

    for ax, (idx, row) in zip(axes, loci_df.iterrows()):
        pts = depth_data[idx]
        if pts:
            x_vals, y_vals = zip(*pts)
            ax.plot(x_vals, y_vals, linewidth=0.7)
        else:
            ax.text(0.5, 0.5, "no coverage", ha="center", va="center", transform=ax.transAxes)

        ax.set_xlim(row.start, row.end)
        ax.ticklabel_format(style="plain", axis="x")
        ax.set_title(row["name"], fontsize=9, pad=2)
        ax.set_xlabel(f"{row.chrom}:{row.start}-{row.end}")
        ax.set_ylabel("depth")

    for ax in axes:
        ax.set_ylim(0, 500)
        ax.set_yticks([0, 250, 500])

    plt.suptitle(depth_gz.name, fontsize=11)
    plt.savefig(out_png, dpi=300)
    plt.close(fig)


# ----- cumulative coverage plot -----

def depths_per_locus(depth_gz: Path, loci_df: pd.DataFrame) -> Iterable[Tuple[str, List[int]]]:
    by_chr: dict[str, List[Tuple[int, int, int]]] = {}
    for idx, r in loci_df.iterrows():
        by_chr.setdefault(r.chrom, []).append((idx, r.start, r.end))
    depth_lists: List[List[int]] = [[] for _ in range(len(loci_df))]

    with gzip.open(depth_gz, "rt") as fh:
        for ln in fh:
            chrom, s, e, d = ln.strip().split("\t")
            s = int(s)
            e = int(e)
            d = int(d)
            for idx, L, R in by_chr.get(chrom, []):
                if e > L and s < R:
                    depth_lists[idx].extend([d] * (min(e, R) - max(s, L)))

    for idx, lst in enumerate(depth_lists):
        yield loci_df.loc[idx, "name"], lst


def cdf(depths: List[int]) -> Tuple[List[int], List[float]]:
    if not depths:
        return [], []
    cnt: dict[int, int] = {}
    total = len(depths)
    rem = total
    for d in depths:
        cnt[d] = cnt.get(d, 0) + 1
    x: List[int] = []
    y: List[float] = []
    for d in sorted(cnt):
        x.append(d)
        y.append(100.0 * rem / total)
        rem -= cnt[d]
    return x, y


def plot_cdf(depth_gz: Path, loci_df: pd.DataFrame, out_png: Path) -> None:
    plt.figure(figsize=(5, 4))
    for name, depths in depths_per_locus(depth_gz, loci_df):
        x, y = cdf(depths)
        if x:
            plt.plot(x, y, label=name, lw=1)

    ax = plt.gca()
    ax.set_xscale("log")
    sf = ScalarFormatter()
    sf.set_scientific(False)
    sf.set_powerlimits((0, 0))
    ax.xaxis.set_major_formatter(sf)
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_yticklabels(["0", "25", "50", "75", "100"])
    ax.set_xlabel("HiFi Coverage (log)")
    ax.set_ylabel("% of the Locus")
    plt.legend(fontsize="xx-small", bbox_to_anchor=(1.02, 0.5), loc="center left", borderaxespad=0)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


# ----- entry point for CLI -----

def plotcov(depth_gz: str, loci: str, out_prefix: str) -> None:
    """Generate facet and cumulative coverage plots."""
    loci_df = load_loci(Path(loci))
    depth_data = collect_depth(Path(depth_gz), loci_df)
    plot_depth_facets(Path(depth_gz), loci_df, depth_data, Path(f"{out_prefix}_bedgraph.png"))
    plot_cdf(Path(depth_gz), loci_df, Path(f"{out_prefix}_covgraph.png"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--depth", required=True)
    p.add_argument("--loci", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    plotcov(args.depth, args.loci, args.out)
