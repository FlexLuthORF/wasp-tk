Installation:

git clone --branch container https://github.com/Watson-IG/wasp.git
cd wasp
bash setup.sh --ref

note: if you want to provide your own reference and bed files, then omit '--ref' option. paths an be set in config. default is /data


Usage:

mapping workflow expects a fofn with sampleId in column 1 and path_to_ccs_pbbam in column 2.
	example:
	nextflow run wasp/main.nf -profile standard --fofn {fofn_path}