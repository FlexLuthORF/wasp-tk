#!/bin/bash
set -e -x

user=$(whoami)
input_file=$(ls *fofn* | head -n 1)
cat $input_file | while read sample ccs
do
    outdir=$PWD/run_hifiasm/${sample}
    mkdir -p $PWD/run_hifiasm/${sample}
    #REMOVE PATH
    sbatch --time=88:00:00 -p compute -o ${outdir}/job.txt --wrap="sh /home/zmvanw01/projects/EF/240505/pipeline.sh ${outdir} fake 10"
    count=`squeue | grep $user | wc -l`
    while [ ${count} -gt 12 ]
    do
	sleep 20s
	count=`squeue | grep $user | wc -l`
    done
done
