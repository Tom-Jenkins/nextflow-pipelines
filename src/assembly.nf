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






