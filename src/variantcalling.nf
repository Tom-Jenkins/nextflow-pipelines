#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.ref_genome = "${PWD}"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to the console
log.info """\
         V A R I A N T - C A L L I N G - N F   P I P E L I N E
         ===================================
         Input reads directory: ${params.reads}
         Input reference genome: ${params.ref_genome}
         Output directory: ${params.outdir}
         Number of threads: ${params.cpus}
         """
         .stripIndent()


// Define workflow
workflow {

    // Import reads
    reads_ch = Channel
        .fromFilePairs("${params.reads}/*_trim_R{1,2}.fq.gz", checkIfExists: false)
        .ifEmpty { error "No paired reads matching the pattern `*_trim_R{1,2}.fq.gz`"}
        .view()

    // Align reads to reference genome


    // See Bioinformatics Workbook Example for intermediate steps

    // Call variants using Freebayes
}

