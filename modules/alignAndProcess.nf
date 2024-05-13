// modules/alignAndProcess.nf

process alignAndProcess {
    publishDir "${params.outdir}/${sampleId}/alg_asm20_to_ref_with_secondarySeq", mode: 'copy'
    //container 'hifi_container.sif'
    
    input:
    tuple val(sampleId), path(merged_reads)

    output:
    tuple val(sampleId), path("${sampleId}.sorted.bam"), emit: final_bam
    tuple val(sampleId), path("${sampleId}.sorted.bam.bai"), emit: final_bam_index

    script:
    """
    minimap2 -x asm20 --secondary-seq -t ${params.cpusPerNode} -L -a ${params.reffn} ${merged_reads} > ${sampleId}.sam
    samtools view -Sbh ${sampleId}.sam > ${sampleId}.bam
    samtools sort -@ ${params.cpusPerNode} ${sampleId}.bam -o ${sampleId}.sorted.bam
    samtools index ${sampleId}.sorted.bam
    """
}
