// modules/extractReads.nf

process extractReads {
    publishDir "${params.outdir}/${sampleId}", mode: 'copy'
    container 'hifi_container.sif'
    
    input:
    tuple val(sampleId), path(ccsPath)

    output:
    tuple val(sampleId), path("${sampleId}.fasta"), emit: reads_fasta
    tuple val(sampleId), path("${sampleId}.fasta.fai"), emit: reads_fasta_fai

    script:
    """
    samtools view ${ccsPath} | awk '{ print ">"$1"\\n"$10 }' > ${sampleId}.fasta
    samtools faidx ${sampleId}.fasta
    """
}