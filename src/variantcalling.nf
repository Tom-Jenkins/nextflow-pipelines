#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.illumina_reads = "${PWD}"
params.ref_genome = "${PWD}"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to the console
log.info """\
         V A R I A N T - C A L L I N G - N F   P I P E L I N E
         ===================================
         Input Illumina directory: ${params.illumina_reads}
         Input reference genome: ${params.nano_reads}
         Output directory: ${params.outdir}
         Number of threads: ${params.cpus}
         """
         .stripIndent()


// Define workflow


