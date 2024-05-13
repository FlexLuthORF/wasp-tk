// modules/alignContigs.nf

process alignContigs {
    publishDir "${params.outdir}/${sampleId}/hifiasm", mode: 'copy'
    container 'hifi_container.sif'
    
    input:
    tuple val(sampleId), path(hap_contigs), path(hap_contigs_fai)

    output:
    tuple val(sampleId), path("asm.bp.hap*.p_ctg_to_ref.sorted.bam"), emit: aligned_contigs

    script:
    """
    for i in 1 2; do
        minimap2 -x asm20 -t ${params.cpusPerNode} -L -a ${params.reffn} asm.bp.hap\${i}.p_ctg.fasta > asm.bp.hap\${i}.p_ctg_to_ref.sam
        samtools view -Sbh asm.bp.hap\${i}.p_ctg_to_ref.sam > asm.bp.hap\${i}.p_ctg_to_ref.bam
        samtools sort -@ ${params.cpusPerNode} asm.bp.hap\${i}.p_ctg_to_ref.bam -o asm.bp.hap\${i}.p_ctg_to_ref.sorted.bam
        samtools index asm.bp.hap\${i}.p_ctg_to_ref.sorted.bam
    done
    """
}