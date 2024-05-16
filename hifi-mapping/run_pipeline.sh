#!/bin/bash
set -e -x

user=$(whoami)
input_file=$1
cat $input_file | while read sample ccs
do
    outdir=$PWD/run_hifiasm/${sample}
    mkdir -p $PWD/run_hifiasm/${sample}
    #REMOVE PATH
    sbatch --time=88:00:00 -p compute -o ${outdir}/job.txt --wrap="sh /home/zmvanw01/git_repos/wasp/hifi-mapping/pipeline.sh ${outdir} ${ccs} 12 ${sample}"
    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 12 ]
    do
	sleep 20s
	count=`squeue | grep $user | wc -l`
    done
done
