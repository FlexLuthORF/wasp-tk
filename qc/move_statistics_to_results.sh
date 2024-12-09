#!/bin/bash

sample=$1
orig_outdir=$2
outdir=$PWD/results/${sample}
statistics_dir=$PWD/statistics/${sample}

mv $statistics_dir/ccs_hifi_readlength_hist_plot.png ${outdir}/stats/read_length_histogram.png
ln -s ${outdir}/stats/read_length_histogram.png $statistics_dir/ccs_hifi_readlength_hist_plot.png
