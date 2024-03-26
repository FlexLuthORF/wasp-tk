import os
import tempfile
import sys


def process_sample(sample_id, r1_gz_path, outdir):
    r2_gz_path = r1_gz_path.replace('R1', 'R2')

    if not os.path.exists(r2_gz_path):
        print(f"Missing R2 file for sample ID {sample_id}")
        return

    with tempfile.NamedTemporaryFile(mode='w+b', suffix='.fastq', delete=False) as r1_tempfile:
        os.system(f"gunzip -c {r1_gz_path} > {r1_tempfile.name}")
        r1_path = r1_tempfile.name

    with tempfile.NamedTemporaryFile(mode='w+b', suffix='.fastq', delete=False) as r2_tempfile:
        os.system(f"gunzip -c {r2_gz_path} > {r2_tempfile.name}")
        r2_path = r2_tempfile.name

    sample_output_dir = os.path.join(outdir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    log_dir = os.path.join(sample_output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    fastqc_dir = os.path.join(log_dir, "fastqc")
    os.makedirs(fastqc_dir, exist_ok=True)
    #if not os.path.exists(f"{sample_output_dir}/CRR_quality-pass.fastq"):
    os.system(f"fastqc -o {log_dir}/fastqc {r2_path}")
    os.system(f"fastqc -o {log_dir}/fastqc {r1_path}")
    if not os.path.exists(f"{sample_output_dir}/CRR_quality-pass.fastq"):      
        os.system(f"FilterSeq.py quality -s {r1_path} -q 20 --nproc 11 --outname CRR --outdir {sample_output_dir} --log {log_dir}/quality-crr.log")

    if not os.path.exists(f"{sample_output_dir}/VRR_quality-pass.fastq"):        
        os.system(f"FilterSeq.py quality -s {r2_path} -q 20 --nproc 11 --outname VRR --outdir {sample_output_dir} --log {log_dir}/quality-vrr.log")

    if os.path.exists(r1_path):   
        os.remove(r1_path)

    if os.path.exists(r2_path):    
        os.remove(r2_path)

    if not os.path.exists(f"{sample_output_dir}/VRR_primers-pass.fastq"):     
        os.system(f"MaskPrimers.py extract -s {sample_output_dir}/VRR_quality-pass.fastq --start 12 --len 4 --barcode --bf BARCODE --mode cut --log {log_dir}/primers-vrr.log --outname VRR --outdir {sample_output_dir} --nproc 11")
    if not os.path.exists(f"{sample_output_dir}/CRR_primers-pass.fastq"):    
          #needs full path to primers
        os.system(f"MaskPrimers.py align -s {sample_output_dir}/CRR_quality-pass.fastq -p /home/zmvanw01/data/references/Human_IG_CRegion_RC.fasta --maxlen 100 --maxerror 0.3 --mode cut --skiprc --pf C_CALL --log {log_dir}/cregion.log --outname CRR --nproc 11")
    if not os.path.exists(f"{sample_output_dir}/table.tab"):    
        os.system(f"ParseLog.py -l {log_dir}/cregion.log -f ID PRIMER ERROR PRSTART --outdir {log_dir}")
    if not os.path.exists(f"{sample_output_dir}/CRR_primers-pass_pair-pass.fastq"):    
        os.system(f"PairSeq.py -1 {sample_output_dir}/VRR_primers-pass.fastq -2 {sample_output_dir}/CRR_primers-pass.fastq --1f BARCODE --2f C_CALL --coord illumina")
    if not os.path.exists(f"{sample_output_dir}/VRR_consensus-pass.fastq"):
        os.system(f"BuildConsensus.py -s {sample_output_dir}/VRR_primers-pass_pair-pass.fastq --bf BARCODE --pf C_CALL --prcons 0.6 -n 1 -q 0 --maxerror 0.1 --maxgap 0.5 --nproc 11 --log {log_dir}/consensus-vrr.log --outdir {sample_output_dir} --outname VRR")
    if not os.path.exists(f"{sample_output_dir}/CRR_consensus-pass.fastq"):   
        os.system(f"BuildConsensus.py -s {sample_output_dir}/CRR_primers-pass_pair-pass.fastq --bf BARCODE --pf C_CALL --prcons 0.6 -n 1 -q 0 --maxerror 0.1 --maxgap 0.5 --nproc 11 --log {log_dir}/consensus-crr.log --outdir {sample_output_dir} --outname CRR")
    if not os.path.exists(f"{sample_output_dir}/CRR_consensus-pass_pair-pass.fastq"):  
        os.system(f"PairSeq.py -1 {sample_output_dir}/VRR_consensus-pass.fastq -2 {sample_output_dir}/CRR_consensus-pass.fastq --coord presto")
    if not os.path.exists(f"{sample_output_dir}/S5_assemble-pass.fastq"):   
        os.system(f"AssemblePairs.py align -1 {sample_output_dir}/VRR_consensus-pass_pair-pass.fastq -2 {sample_output_dir}/CRR_consensus-pass_pair-pass.fastq --coord presto --rc tail --1f CONSCOUNT --2f PRCONS CONSCOUNT --minlen 8 --maxerror 0.3 --alpha 1e-5 --scanrev --nproc 11 --log {log_dir}/assemble.log --outname S5")
    if not os.path.exists(f"{sample_output_dir}/S5-MQ_maskqual-pass.fastq"):    
        os.system(f"FilterSeq.py maskqual -s {sample_output_dir}/S5_assemble-pass.fastq -q 30 --nproc 11 --outname S5-MQ --log {log_dir}/maskqual.log")
    if not os.path.exists(f"{sample_output_dir}/S5-final_reheader.fastq"):    
        os.system(f"ParseHeaders.py collapse -s {sample_output_dir}/S5-MQ_maskqual-pass.fastq -f CONSCOUNT --act min --outname S5-final")
    if not os.path.exists(f"{sample_output_dir}/S5-final_total.fastq"):    
        os.system(f"mv {sample_output_dir}/S5-final_reheader.fastq {sample_output_dir}/S5-final_total.fastq")
    if not os.path.exists(f"{sample_output_dir}/S5-final_total_collapse-unique.fastq"):    
        os.system(f"CollapseSeq.py -s {sample_output_dir}/S5-final_total.fastq -n 0 --inner --uf PRCONS --cf CONSCOUNT --act sum ")
    if not os.path.exists(f"{sample_output_dir}/S5-final_total_collapse-unique_atleast-2.fasta"):    
        os.system(f"SplitSeq.py group -s {sample_output_dir}/S5-final_total_collapse-unique.fastq -f CONSCOUNT --num 2")


    os.system(f'ParseHeaders.py table -s "{sample_output_dir}/S5-final_total.fastq" -f ID PRCONS CONSCOUNT --outname "final-total" --outdir {sample_output_dir}/logs')
    os.system(f'ParseHeaders.py table -s "{sample_output_dir}/S5-final_total_collapse-unique.fastq" -f ID PRCONS CONSCOUNT DUPCOUNT --outname "final-unique" --outdir {sample_output_dir}/logs')
    os.system(f'ParseHeaders.py table -s "{sample_output_dir}/S5-final_total_collapse-unique_atleast-2.fastq" -f ID PRCONS CONSCOUNT DUPCOUNT --outname "final-unique-atleast2" --outdir {sample_output_dir}/logs')

    os.system(f'ParseLog.py -l "{sample_output_dir}/logs/primers-vrr.log" -f ID BARCODE ERROR --outdir {sample_output_dir}/logs')
    os.system(f'ParseLog.py -l "{sample_output_dir}/logs/consensus-vrr.log" "{sample_output_dir}/logs/consensus-crr.log" -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR --outdir {sample_output_dir}/logs')
    os.system(f'ParseLog.py -l "{sample_output_dir}/logs/assemble.log" -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 --outdir {sample_output_dir}/logs')
    os.system(f'ParseLog.py -l "{sample_output_dir}/logs/maskqual.log" -f ID MASKED --outdir {sample_output_dir}/logs')

# changeo starts
    if not os.path.exists(f"{sample_output_dir}/S5-final_total_collapse-unique_atleast-2_reheader.fasta"):
        os.system(f'ParseHeaders.py rename -s {sample_output_dir}/S5-final_total_collapse-unique_atleast-2.fastq --fasta -f PRCONS -k C_CALL')
    if not os.path.exists(f"{sample_output_dir}/S5_igblast.fmt7"):
        os.system(f'AssignGenes.py igblast -s {sample_output_dir}/S5-final_total_collapse-unique_atleast-2_reheader.fasta --organism human --loci ig -b /home/zmvanw01/share/igblast --format blast --nproc 11 --outdir {sample_output_dir} --outname "S5"')
    if not os.path.exists(f"{sample_output_dir}/S5_db-pass.tsv"):
        os.system(f'MakeDb.py igblast -s {sample_output_dir}/S5-final_total_collapse-unique_atleast-2_reheader.fasta -i {sample_output_dir}/S5_igblast.fmt7 --extended --failed --format airr -r /home/zmvanw01/share/germlines/imgt/human/vdj/ --outname S5')

if __name__ == '__main__':
    # Keeping the first components up to the required depth and appending 'presto/{sample_id}'
    sample_id = sys.argv[1]
    r1_gz_path = sys.argv[2]
    base_dir = os.path.normpath(r1_gz_path).split(os.sep)
    outdir = "./presto"
    process_sample(sample_id, r1_gz_path, outdir)