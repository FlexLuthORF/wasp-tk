#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Process imports
include { extractReads } from './modules/extractReads'
include { runHifiasm } from './modules/runHifiasm'
include { alignContigs } from './modules/alignContigs'
include { alignUnitigs } from './modules/alignUnitigs'
include { extractSoftClip } from './modules/extractSoftClip'
include { mergeAndDedup } from './modules/mergeAndDedup'
include { alignAndProcess } from './modules/alignAndProcess'

// Main workflow
workflow {
    // Read the FOFN file and create input channel
    ccsData = Channel.fromPath(params.fofn)
                     .splitCsv(sep: '\t')
                     .map { sampleId, ccsPath -> [sampleId, file(ccsPath)] }

    // Run the extractReads process
    extractReads(ccsData)

    // Run the runHifiasm process
    runHifiasm(extractReads.out.reads_fasta, extractReads.out.reads_fasta_fai)

    // Run the alignContigs process
    alignContigs(runHifiasm.out.hap_contigs)

    // Run the alignUnitigs process
    alignUnitigs(runHifiasm.out.unitigs)

    // Run the extractSoftClip process
    extractSoftClip(alignContigs.out.aligned_contigs)

    // Run the mergeAndDedup process
    mergeAndDedup(extractSoftClip.out.softclip_aligned)

    // Run the alignAndProcess process
    alignAndProcess(mergeAndDedup.out.merged_reads)
}