// modules/mergeAndDedup.nf

process mergeAndDedup {
    publishDir "${params.outdir}/${sampleId}/merged_bam", mode: 'copy'
    container 'hifi_container.sif'
    
    input:
    tuple val(sampleId), path(aligned_contigs)

    output:
    tuple val(sampleId), path("merged_all_reads.rmdup.fasta"), emit: merged_reads

    script:
    """
    samtools merge -f merged.bam ${aligned_contigs}
    samtools sort -@ ${params.cpusPerNode} merged.bam -o merged.sorted.bam
    samtools index merged.sorted.bam
    samtools fasta --reference ${params.reffn} merged.sorted.bam > merged_all_reads.fasta
    seqkit rmdup --by-seq merged_all_reads.fasta -o merged_all_reads.rmdup.fasta
    """
}