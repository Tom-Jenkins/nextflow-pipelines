#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.readsPattern = "${params.reads}/*_R{1,2}_001.fastq.gz"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to console
println "Input directory: ${params.reads}"
println "Output directory: ${params.outdir}/trimmed_reads"
println "Number of threads: ${params.cpus}"

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
    // Use replaceAll() to:
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
    reads_ch = Channel.fromFilePairs(params.readsPattern, checkIfExists: true)
    FASTP(reads_ch)
}
