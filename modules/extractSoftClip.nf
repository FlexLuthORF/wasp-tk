// modules/extractSoftClip.nf

process extractSoftClip {
    container 'hifi_container.sif'
    publishDir "${params.outdir}/${sampleId}/break_at_soft_clip", mode: 'copy'

    input:
    tuple val(sampleId), path(aligned_contigs)

    output:
    tuple val(sampleId), path("*_hifi_asm.fasta"), path("*_hifi_asm.fasta.fai"), emit: softclip_fasta
    tuple val(sampleId), path("*_hifi_asm_to_ref.sorted.bam"), path("*_hifi_asm_to_ref.sorted.bam.bai"), emit: softclip_aligned

    script:
    """
    for i in 1 2; do
        extract_soft_clip_seq.py asm.bp.hap\${i}.p_ctg_to_ref.sorted.bam > \${i}_hifi_asm.fasta
        samtools faidx \${i}_hifi_asm.fasta
        
        minimap2 -t ${params.cpusPerNode} -L -a ${params.reffn} \${i}_hifi_asm.fasta > \${i}_hifi_asm_to_ref.sam
        samtools view -Sbh \${i}_hifi_asm_to_ref.sam > \${i}_hifi_asm_to_ref.bam
        samtools sort -@ ${params.cpusPerNode} \${i}_hifi_asm_to_ref.bam -o \${i}_hifi_asm_to_ref.sorted.bam
        samtools index \${i}_hifi_asm_to_ref.sorted.bam
        
        minimap2 -x asm20 -t ${params.cpusPerNode} -L -a ${params.reffn} \${i}_hifi_asm.fasta > \${i}_asm20_hifi_asm_to_ref.sam
        samtools view -Sbh \${i}_asm20_hifi_asm_to_ref.sam > \${i}_asm20_hifi_asm_to_ref.bam
        samtools sort -@ ${params.cpusPerNode} \${i}_asm20_hifi_asm_to_ref.bam -o \${i}_asm20_hifi_asm_to_ref.sorted.bam
        samtools index \${i}_asm20_hifi_asm_to_ref.sorted.bam
    done
    """
}