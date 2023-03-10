#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.sra = "${PWD}"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to the console
log.info """\
         F A S T P - N F   P I P E L I N E
         ===================================
         Input directory: ${params.reads}
         Output directory: ${params.outdir}/trimmed_reads
         Number of threads: ${params.cpus}
         """
         .stripIndent()

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
    // 1. Remove "10628_" (project ID) from file name
    // 2. Remove "S[1 or more number(s)]" (plate ID) from file name
    tuple val(sample_id), path("${sample_id.replaceAll(/10628_/, "").replaceAll(/_S[0-9]*/, "")}_trim_R*.fq.gz")

    script:
    """
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${sample_id.replaceAll(/10628_/, "").replaceAll(/_S[0-9]*/, "")}_trim_R1.fq.gz \
    -O ${sample_id.replaceAll(/10628_/, "").replaceAll(/_S[0-9]*/, "")}_trim_R2.fq.gz \
    --cut_tail \
    --cut_tail_window_size 4 \
    --cut_tail_mean_quality 22 \
    --detect_adapter_for_pe \
    --trim_poly_g \
    --length_required 75 \
    --dont_eval_duplication \
    --thread ${task.cpus}
    """
}

// Define workflow
workflow {

    reads_ch = Channel
        .fromFilePairs(["${params.reads}/*_R{1,2}_001.fastq.gz", "${params.sra}/*_{1,2}.fastq.gz"], checkIfExists: false)
        .ifEmpty { error "No paired reads matching the pattern `*_R{1,2}_001.fastq.gz` or `*_{1,2}.fastq.gz` at ${params.reads}" }
        // .view()

    FASTP(reads_ch)
}
