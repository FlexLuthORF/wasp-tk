#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def currentDate = new Date().format('yyyyMMdd')

params.outdir = "${System.getProperty('user.dir')}/${currentDate}_output"
params.fofn = params.fofn
params.cpusPerNode = params.cpusPerNode
params.reffn = params.reffn
params.maxJobs = params.maxJobs

include { extractReads } from './modules/extractReads'
include { runHifiasm } from './modules/runHifiasm'
include { alignContigs } from './modules/alignContigs'
//include { alignUnitigs } from './modules/alignUnitigs'
include { extractSoftClip } from './modules/extractSoftClip'
include { mergeAndDedup } from './modules/mergeAndDedup'
include { alignAndProcess } from './modules/alignAndProcess'

workflow {
    ccsData = Channel.fromPath(params.fofn)
                     .splitCsv(sep: '\t', header: false)
                     .map { sample_id, ccs_path -> [sample_id, file(ccs_path)] }

    extractReads(ccsData)
    runHifiasm(extractReads.out.reads_fasta)
    
    // Collect outputs from runHifiasm
    hap_contigs_hap1 = runHifiasm.out.hap_contigs_hap1
    hap_contigs_hap2 = runHifiasm.out.hap_contigs_hap2
    
    // Pass the collected outputs to alignContigs
    alignContigs(hap_contigs_hap1, hap_contigs_hap2)
    extractSoftClip(alignContigs.out.aligned_contigs)
    mergeAndDedup(extractSoftClip.out.softclip_aligned)
    alignAndProcess(mergeAndDedup.out.merged_reads)
}