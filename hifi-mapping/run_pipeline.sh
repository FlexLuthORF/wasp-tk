#!/bin/bash
set -e -x

user=$(whoami)
cat \*fofn\* | while read sample ccs
do
    outdir=$PWD/run_hifiasm/${sample}
    mkdir -p $PWD/run_hifiasm/${sample}
    sbatch --time=88:00:00 -p compute -o ${outdir}/job.txt --wrap="sh pipeline.sh ${outdir} ${ccs} 10"
    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 18 ]
    do
	sleep 1s
	count=`squeue | grep $user | wc -l`
    done
done
