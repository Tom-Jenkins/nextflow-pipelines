#!/usr/bin/env nextflow

// Parameters
params.reads = "${PWD}"
params.suffix = "_{1,2}.fastq.gz"
params.transcriptome = "${PWD}"
params.salmon = false
params.salmonParams = "--libType A --gcBias --seqBias"
params.kallisto = false
params.kallistoParams = " "
// params.star = false
// params.starIndex = "${PWD}"
// params.starParams = ""
params.outdir = "${PWD}"
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         Q U A N T I F Y E X P R E S S I O N - N F   S C R I P T
         ===================================
         Input reads directory: ${params.reads}
         Input reads suffix: ${params.suffix}
         Transcriptome FASTA: ${params.transcriptome}
         Salmon parameter(s): ${params.salmon ? params.salmonParams : "n/a"}
         Kallisto parameter(s): ${params.kallisto ? params.kallistoParams : "n/a"}
         Output directory: ${params.outdir}
         Number of threads: ${params.cpus}
         Script version: v0.2
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

        // Salmon
        if (params.salmon) {

            // Index reference
            salmon_index = SALMON_INDEX(params.transcriptome)

            // Quantify expression
            SALMON_QUANT(salmon_index, reads_ch)
        }

        // Kallisto
        if (params.kallisto) {

            // Index reference
            kallisto_index = KALLISTO_INDEX(params.transcriptome)

            // Quantify expression
            KALLISTO_QUANT(kallisto_index, reads_ch)
        }

        // STAR
        // if (params.star) {

        //     // Align reads
        //     STAR_ALIGN(reads_ch)
        // }
    }
}


// Index reference for Salmon
process SALMON_INDEX {

    // Directives
    maxForks 1 // set maximum number of parallel tasks to 1

    input:
    path(transcriptome)

    output:
    path("index")

    script:
    """
    salmon index \
        -t ${transcriptome} \
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

// Index reference for Kallisto
process KALLISTO_INDEX {

    // Directives
    maxForks 1 // set maximum number of parallel tasks to 1

    input:
    path(transcriptome)

    output:
    path("index")

    script:
    """
    kallisto index \
        --index=index \
        --threads=${params.cpus} \
        ${transcriptome}
    """
}

// Quantify expression using Kallisto
process KALLISTO_QUANT {

    // Directives
    errorStrategy "ignore"
    maxForks 1
    publishDir "${params.outdir}/kallisto_output", mode: "copy"

    input:
    path(index)
    tuple val(sample_id), path(reads)

    output:
    path("*")

    script:
    """
    kallisto quant \
        --index=${index} \
        ${params.kallistoParams} \
        --output-dir=${sample_id} \
        --threads=${params.cpus} \
        ${reads[0]} ${reads[1]}
    """
}

// Align reads using STAR (testing)
process STAR_ALIGN {

    // Directives
    maxForks 1 // set maximum number of parallel tasks to 1
    publishDir "${params.outdir}/STAR_output", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}.bam")
    path("${sample_id}.counts"), optional: true

    script:
    """
    STAR \
        --runMode alignReads \
        --genomeDir ${params.starIndex} \
        ${params.starParams} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --outFileNamePrefix ${sample_id}_ \
        --runThreadN ${params.cpus}
    
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    mv ${sample_id}_ReadsPerGene.out.tab ${sample_id}.counts     
    """
}