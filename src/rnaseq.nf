#!/usr/bin/env nextflow

// Parameters
params.reads = "${PWD}"
params.suffix = "_{1,2}.fastq.gz"
params.reference = "${PWD}"
params.annotation = "${PWD}"
params.outdir = "${PWD}"
params.salmon = false
params.salmonParams = "--libType A --gcBias --seqBias"
params.star = false
params.starParams = ""
params.featureCounts = false
params.featureCountsParams = ""
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         R N A S E Q - N F   S C R I P T
         ===================================
         Input reads directory: ${params.reads}
         Input reads suffix: ${params.suffix}
         Reference FASTA: ${params.reference}
         Annotation GTF: ${params.annotation}
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

        // STAR
        if (params.star) {

            // Index reference
            star_index = STAR_INDEX(params.reference, params.annotation)

            // Align reads
            // STAR_ALIGN(star_index, reads_ch)
        }

        // Salmon
        if (params.salmon) {

            // Index reference
            salmon_index = SALMON_INDEX(params.reference)

            // Quantify expression
            SALMON_QUANT(salmon_index, reads_ch)
        }
    }
}


// Index reference for STAR
process STAR_INDEX {

    // Directives
    maxForks 1 // set maximum number of parallel tasks to 1

    input:
    path(reference)
    path(annotation)

    output:
    path("index")

    script:
    """
    STAR \
        --runMode genomeGenerate \
        --genomeFastaFiles ${reference} \
        --sjdbGTFfile ${annotation} \
        --runThreadN ${params.cpus}
    """
}

// Index reference for Salmon
process SALMON_INDEX {

    // Directives
    maxForks 1 // set maximum number of parallel tasks to 1

    input:
    path(reference)

    output:
    path("index")

    script:
    """
    salmon index \
        -t ${reference} \
        -i index \
        --threads ${params.cpus}
    """
}

// Quantify expression using Salmon
process SALMON_QUANT {

    // Directives
    errorStrategy "ignore"
    maxForks 1
    publishDir "${params.outdir}/salmon_output", mode: "copy"

    input:
    path(index)
    tuple val(sample_id), path(reads)

    output:
    path("*")

    script:
    """
    salmon quant \
        -i ${index} \
        ${params.salmonParams} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        --output ${sample_id} \
        --threads ${params.cpus}
    """
}
