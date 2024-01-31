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
        // .view()

    // Align reads to reference genome
    bam_ch = ALIGN_TO_REF_GENOME(reads_ch)

    // Process bam files and pipe the output to bcftools v1.17
    bam_ch | PROCESS_BAM | collect | CALL_VARIANTS
    
}


process ALIGN_TO_REF_GENOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"

    // [sample_ID, [sample_ID_trim.R1.fq.gz, sample_ID_trim.R2.fq.gz]]
    input:
    tuple val(sample_id), path(reads)

    // [sample_ID, [sample_ID.bam]]
    output:
    tuple val(sample_id), path("${sample_id}.bam"), optional: true

    // Index reference genome
    // Align reads using bwa-mem2
    // Filter out unmapped reads and convert to bam
    script:
    """
    bwa-mem2 index ${params.ref_genome}
    bwa-mem2 mem -M -t ${params.cpus} -a ${params.ref_genome} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    samtools view -@ ${params.cpus} -F 4 -b ${sample_id}.sam > ${sample_id}.bam
    rm ${sample_id}.sam
    """
}

process PROCESS_BAM {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    publishDir "${params.outdir}/processed_bams", mode: "copy"

    // [sample_ID, [sample_ID.bam]]
    input:
    tuple val(sample_id), path(bam)

    output:
    path("*-sorted-md-rg.bam")

    // Sort bam
    // Mark duplicates
    // Add read groups
    // Index bam
    script:
    """
    java -jar ~/software/picard.jar SortSam \
        I=${bam} \
        O=${sample_id}-sorted.bam \
        SORT_ORDER=coordinate

    java -jar ~/software/picard.jar MarkDuplicates \
        I=${sample_id}-sorted.bam \
        O=${sample_id}-sorted-md.bam \
        M=${sample_id}-md-metrics.txt

    java -jar ~/software/picard.jar AddOrReplaceReadGroups \
        I=${sample_id}-sorted-md.bam \
        O=${sample_id}-sorted-md-rg.bam \
        RGID=${sample_id} \
        RGLB=${sample_id} \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=${sample_id}

    rm ${sample_id}-sorted.bam ${sample_id}-sorted-md.bam
    samtools index ${sample_id}-sorted-md-rg.bam
    """
}

process CALL_VARIANTS {

    // Directives
    cpus params.cpus
    publishDir "${params.outdir}", mode: "copy"

    input:
    path(bam)

    output:
    path("bcftools.vcf")

    // Freebayes v1.3.6 ran for too long (even in parallel mode)
    // freebayes \
    //     --fasta-reference ${params.ref_genome} \
    //     --bam ${bam} \
    //     --ploidy 2 \
    //     --vcf freebayes_unfiltered.vcf

    // Run bcftools mpileup
    // Run bcftools call
    script:
    """
    bcftools mpileup \
        --threads ${params.cpus} \
        --fasta-ref ${params.ref_genome} \
        --min-MQ 20 \
        --min-BQ 30 \
        --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
        --output bcftools.pileup \
        --output-type z \
        ${bam}
    
    bcftools call \
        --threads ${params.cpus} \
        --multiallelic-caller \
        --variants-only \
        --ploidy 2 \
        --output-type v \
        --output bcftools.vcf \
        bcftools.pileup
    """
}
