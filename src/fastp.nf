#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.sra = "${PWD}"
params.outdir = "${PWD}"
params.project_code = ""
params.cpus = 16

// Print parameters to the console
log.info """\
         F A S T P - N F   P I P E L I N E
         ===================================
         Input directory: ${params.reads}
         Output directory: ${params.outdir}/trimmed_reads
         Project code(s): ${params.project_code}
         Number of threads: ${params.cpus}
         """
         .stripIndent()

// Define workflow
workflow {

    // Import reads
    reads_ch = Channel
        .fromFilePairs(["${params.reads}/*_R{1,2}_001.fastq.gz", "${params.sra}/*_{1,2}.fastq.gz"], checkIfExists: false)
        .ifEmpty { error "No paired reads matching the pattern `*_R{1,2}_001.fastq.gz` or `*_{1,2}.fastq.gz` at ${params.reads}" }
        // .view()

    // Trim reads
    FASTP(reads_ch)
}


// FASTP
process FASTP {

    // Directives
    cpus params.cpus
    publishDir "${params.outdir}/trimmed_reads", mode: "copy"
    conda "/lustre/home/tj311/software/mambaforge3/envs/fastp"
    // conda "fastp=0.23.2"

    input:
    tuple val(sample_id), path(reads)

    output:
    // Output file names in the format `sampleID_trim.fq.gz`
    // This code only gets executed for the files provided by the Sequencing Facility (it should not execute on SRA files)
    // 1. Remove project ID (e.g. "10628|11002") from file name
    // 2. Remove any underscores that remain at the start of file name
    // 3. Remove "S[1 or more number(s)]" (plate ID) from file name
    tuple val(sample_id), path("${sample_id.replaceAll(/${params.project_code}_/, "").replaceFirst(/^_/, "").replaceAll(/_S[0-9]*/, "")}_trim_R*.fq.gz")

    script:
    """
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${sample_id.replaceAll(/${params.project_code}_/, "").replaceFirst(/^_/, "").replaceAll(/_S[0-9]*/, "")}_trim_R1.fq.gz \
    -O ${sample_id.replaceAll(/${params.project_code}_/, "").replaceFirst(/^_/, "").replaceAll(/_S[0-9]*/, "")}_trim_R2.fq.gz \
    --qualified_quality_phred 20 \
    --trim_poly_g \
    --length_required 75 \
    --thread ${task.cpus}
    """
}


