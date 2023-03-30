#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.pcalcareumREF = "${PWD}"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to the console
log.info """\
         A S S E M B L Y - N F   P I P E L I N E
         ===================================
         Input directory: ${params.reads}
         Output directory: ${params.outdir}
         Phymatolithon calcareum reference: ${params.pcalcareumREF}
         Number of threads: ${params.cpus}
         """
         .stripIndent()


// Define workflow
workflow {

    // Import Nanopore reads
    reads_ch = Channel
        .fromPath("${params.reads}/*fastq.gz", checkIfExists: false)
        .ifEmpty { error "No reads matching the pattern `*fastq.gz`" }
        .view()

    // Detect and remove adapters using Porechop v0.2.4
    // porechop_ch = PORECHOP(reads_ch)
}

// PORECHOP
process PORECHOP {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz")

    script:
    """
    porechop -i ${reads} -o ${sample_id}/.fq -t ${params.cpus}
    gzip ${sample_id}/.fq
    """
}



