// modules/runHifiasm.nf

process runHifiasm {
    publishDir "${params.outdir}/${sampleId}/hifiasm", mode: 'copy'
    container 'hifi_container.sif'
    
    input:
    tuple val(sampleId), path(reads_fasta), path(reads_fasta_fai)

    output:
    tuple val(sampleId), path("asm.bp.hap*.p_ctg.fasta"), emit: hap_contigs
    tuple val(sampleId), path("asm.bp.hap*.p_ctg.fasta.fai"), emit: hap_contigs_fai

    script:
    """
    hifiasm -o asm -t ${params.cpusPerNode} ${reads_fasta}

    for i in 1 2; do
        gfatools gfa2fa asm.bp.hap\${i}.p_ctg.gfa > asm.bp.hap\${i}.p_ctg.fasta
        samtools faidx asm.bp.hap\${i}.p_ctg.fasta
    done
    """
}