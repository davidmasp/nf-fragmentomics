#!/usr/bin/env nextflow
/*
========================================================================================
FRAGMENTOMICS PIPELINE
========================================================================================
*/

nextflow.enable.dsl = 2

include { FRAGMENTOMICS } from './workflows/fragmentomics.nf'
include { VERSIONS } from './modules/local/utils/versions.nf'
include { FA_TO_TWOBIT } from './modules/local/ucsctools/main.nf'

def create_target_channel(LinkedHashMap row) {
    return [row.source, file(row.bed)]
}

def create_sample_channel(LinkedHashMap row) {
    // create all at once
    def meta = [
        caseid: row.caseid,
        sampleid: row.sampleid,
        timepoint: row.timepoint
    ]

    // bam, bai and bw can be null
    return [
        meta,
        row.bam ? file(row.bam) : null,
        row.bai ? file(row.bai) : null,
        row.bw  ? file(row.bw) : null
    ]
}

// MAIN WORKFLOW
workflow {

    main:
    // Init param files
    genome_fasta = params.genome_fasta ? channel.fromPath(params.genome_fasta) : channel.empty() 
    blacklist_bed = params.blacklist_bed ? channel.fromPath(params.blacklist_bed) : channel.empty()

    ch_versions = channel.empty()
    if (params.genome_fasta) {
        FA_TO_TWOBIT(genome_fasta)
        genome_2bit = FA_TO_TWOBIT.out.twobit
        ch_versions = ch_versions.mix(FA_TO_TWOBIT.out.versions)
    } else {
        genome_2bit = params.genome_2bit ? channel.fromPath(params.genome_2bit) : channel.empty()
    }

    // samples channel
    sample_ch = channel.fromPath(params.input)
        .splitCsv(header:true, sep:',')
        .map{ srow ->  create_sample_channel(srow) }

    // targets channel
    target_ch = channel.fromPath(params.targets)
        .splitCsv(header: true, sep:',')
        .map{ trow -> create_target_channel(trow) }

    // filter targets for lines if not stubrun
    if (workflow.stubRun == false) {
        target_ch = target_ch
            .filter{ it ->
                it[1].readLines().size() > 1
            }
    }


    target_ch = target_ch
        .groupTuple(
            by: 0,
            size: params.collate_size,
            remainder: true)
     
    FRAGMENTOMICS(
        sample_ch,
        target_ch,
        genome_2bit,
        blacklist_bed
    )

    // collect versions in a single file simple mode
    VERSIONS(
        FRAGMENTOMICS.out.versions
            .mix(ch_versions)
            .map { version_file ->
                version_file.text
            }
            .collect()
    )
}
