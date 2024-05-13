// modules/runHifiasm.nf

process runHifiasm {
    publishDir "${params.outdir}/${sampleId}/hifiasm", mode: 'copy'
    //container 'hifi_container.sif'

    input:
    tuple val(sampleId), path(reads_fasta)

    output:
    tuple val(sampleId), path("asm.bp.hap1.p_ctg.fasta"), path("asm.bp.hap1.p_ctg.fasta.fai"), emit: hap_contigs_hap1
    tuple val(sampleId), path("asm.bp.hap2.p_ctg.fasta"), path("asm.bp.hap2.p_ctg.fasta.fai"), emit: hap_contigs_hap2

    script:
    """
    hifiasm -o asm -t ${params.cpusPerNode} ${reads_fasta}

    for i in 1 2; do
        gfatools gfa2fa asm.bp.hap\${i}.p_ctg.gfa > asm.bp.hap\${i}.p_ctg.fasta
        samtools faidx asm.bp.hap\${i}.p_ctg.fasta
    done
    """
}
