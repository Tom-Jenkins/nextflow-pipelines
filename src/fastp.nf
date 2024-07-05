#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Example usage:
// nextflow run ./nextflow-pipelines/src/fastp.nf \
//     --reads /path/to/reads/directory \
//     --suffix "_{1,2}.fastq.gz" \
//     --adapters /path/to/adapters/adapters.fasta \
//     --filter "--qualified_quality_phred 30 --length_required 100 --trim_poly_g" \
//     --outdir /path/to/output/directory \
//     --cpus 16

// Parameters
params.reads = "${PWD}"
params.suffix = "_{1,2}.fastq.gz"
params.adapters = ""
params.outdir = "${PWD}"
params.filter = "--qualified_quality_phred 30 --length_required 100 --trim_poly_g"
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         F A S T P - N F   P I P E L I N E
         ===================================
         Input directory: ${params.reads}
         Output directory: ${params.outdir}/trimmed_reads
         Filtering parameter(s): ${params.filter}
         Number of threads: ${params.cpus}
         Script version: v0.2
         """
         .stripIndent()

// Define workflow
workflow {

    // Create reads channel
    reads_ch = Channel
        .fromFilePairs("${params.reads}/*${params.suffix}", checkIfExists: false)
        .ifEmpty { 
            println "No paired reads found with the pattern `*${params.suffix}`"
        }    

    // Test run to view parameters and contents of reads_ch
    if ( params.test ) {
        reads_ch.view()
    }
    // Filter and trim reads using fastp
    else {
        FASTP(reads_ch)
    }
}

// FASTP
process FASTP {

    // Directives
    cpus params.cpus
    publishDir "${params.outdir}/trimmed_reads", mode: "copy"

    input:
    tuple val(sample_id), path(reads)

    output:
    // Output file names in the format `sample_ID.fp.fq.gz`
    path("*.fp.fq.gz")
    path("*.html")
    // path("*.json")

    script:
    """
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${sample_id}_1.fp.fq.gz \
    -O ${sample_id}_2.fp.fq.gz \
    --adapter_fasta ${params.adapters} \
    ${params.filter} \
    --json ${sample_id}_fastp.json \
    --html ${sample_id}_fastp.html \
    --thread ${task.cpus}
    """
}
