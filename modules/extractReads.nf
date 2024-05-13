// modules/extractReads.nf

process extractReads {
    publishDir "${params.outdir}/${sampleId}", mode: 'copy'
    //container 'hifi_container.sif'

    input:
    tuple val(sampleId), path(ccsPath)

    output:
    tuple val(sampleId), path("${sampleId}.fasta"), emit: reads_fasta

    script:
    """
    samtools view ${ccsPath} | awk '{ print \">\"\$1\"\\n\"\$10 }' > ${sampleId}.fasta
    samtools faidx ${sampleId}.fasta
    """
}
