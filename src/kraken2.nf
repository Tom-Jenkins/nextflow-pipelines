#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.suffix = "_{1,2}.fastq.gz"
params.krakenDB = "${PWD}"
params.kraken2 = ""
params.outdir = "${PWD}"
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         K R A K E N 2 - N F   S C R I P T
         ===================================
         Input reads directory: ${params.reads}
         Input reads suffix: ${params.suffix}
         Kraken2 database: ${params.krakenDB}
         Kraken2 parameter(s): ${params.kraken2}
         Output directory: ${params.outdir}
         Number of threads: ${params.cpus}
         Script version: v0.1
         """
         .stripIndent()


// Define workflow
workflow {

    // Read in paired reads
    reads_ch = Channel
        .fromFilePairs("${params.reads}/*${params.suffix}", checkIfExists: false)
        .ifEmpty { error "No paired reads matching the pattern `sampleID${params.suffix}`"}

    // Test run to view parameters and contents of reads_ch
    if (params.test) {
        reads_ch.view()
    }

    // Run main pipeline
    else {
        // KRAKEN2
        KRAKEN2(reads_ch)
    }
}

// KRAKEN2
process KRAKEN2 {

    // Directives
    publishDir "${params.outdir}", mode: "copy"
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1

    // [sample_ID, [read1, read2]]
    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}.kraken")
    path("${sample_id}.report")

    script:
    """
    kraken2 \
        --threads ${params.cpus} \
        --db ${params.krakenDB} \
        ${params.kraken2} \
        --output ${sample_id}.kraken \
        --report ${sample_id}.report \
        --paired ${reads[0]} ${reads[1]}
    """
}
