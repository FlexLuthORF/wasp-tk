ASSUMES SLURM HPC IN USE

Need to install:

IGenotyper - branch/zach-changes
	https://github.com/oscarlr/IGenotyper/tree/zach-changes

Install minimap2
	https://github.com/lh3/minimap2

Install hifiasm
	https://github.com/chhylp123/hifiasm

Install VDJbase Genomics

Install receptor_utils

Install immcantation env from immcantation.yml


Usage for WGS/capture:

Note: Various steps in the pipeline expect to find *fofn* that has donorID in col1, and path to raw data in col2.

1. User run_pipeline.sh to iterate the mapping and soft clip extraction. *fofn* is expected to be found in directory it is ran in, the reads will be in ./run_hifiasm/{donorID}. This will also merge hap1 and hap2, as well as ensure all contigs have a unique name preparing it for makegenes.py

	bash run_pipeline.sh *fofn*
	
2. Use iter_makegenes.py to iterate the annotations.csv file generatoin using VDJbase-Genomics scripts. accepts *fofn* as an argument

	python iter_makegenes.py *fofn*



Usage for repSeq:

Note: *fofn* used is expected to be donorID col1 and full path to R1 paired read in col2

1. Activate immcantation env

2. Run iterate_presto.py with the *fofn* supplied as an argument. It will create a ./presto folder and fill in donorId folders there

3. Run iter_personal_repseq.py with *fof* supplies as an argument. It expects to be ran in same root as above and will create ./presto/changeo

4. Switch to r enviroment and run defineclones.R with *fofn* supplied as argument. 




Read support usage:

NEEDS ADDED
