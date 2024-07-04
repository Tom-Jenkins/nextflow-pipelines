#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.sampleSheet = "${PWD}/sample_sheet.csv"
params.genome = "${PWD}"
params.outdir = "${PWD}"
params.vcf = "bcftools_variants"
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         V A R I A N T - C A L L I N G - N F   P I P E L I N E
         ===================================
         Input reads sample sheet: ${params.sampleSheet}
         Input reference genome: ${params.genome}
         Output directory: ${params.outdir}
         VCF output prefix: ${params.vcf}
         Number of threads: ${params.cpus}
         """
         .stripIndent()


// Define workflow
workflow {
   
    // Read the CSV file and create a channel of tuples
    sample_ch = Channel
        .fromPath(params.sampleSheet)
        .splitCsv(header: true)
        .map { row ->
            // Extract sample ID, read 1 path, and read 2 path from each row
            def sampleId = row["sample"]
            def library = row["library"]
            def run = row["run"]
            def read1Path = row["read1"]
            def read2Path = row["read2"]

            // Create a tuple with the desired structure
            return [sampleId, library, run, [read1Path, read2Path]]
        }

    // Test run to view parameters and contents of sample_ch
    if ( params.test ) {
        sample_ch.view()
    }
    // Run main pipeline
    else {
        // Align reads to reference genome
        bam_ch = ALIGN_TO_REF_GENOME(sample_ch)

        // Process bam files and pipe the output to bcftools v1.19
        bam_ch | PROCESS_BAM | collect | CALL_VARIANTS
    }    
}


process ALIGN_TO_REF_GENOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"

    input:
    tuple val(sample_id), val(library), val(run), path(reads)

    // [sample_ID, [sample_ID.bam]]
    output:
    tuple val(sample_id), val(library), val(run), path("${sample_id}.sorted.bam"), optional: true

    // Align reads using bowtie2
    // Filter out unmapped reads, convert to bam, and sort bam
    script:
    """
    bowtie2 -p ${task.cpus} -x ${params.genome} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam
    samtools view -F 4 --threads ${task.cpus} -b ${sample_id}.sam > ${sample_id}.bam
    samtools sort --threads ${task.cpus} ${sample_id}.bam > ${sample_id}.sorted.bam
    rm ${sample_id}.sam ${sample_id}.bam
    """
}

process PROCESS_BAM {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    publishDir "${params.outdir}/processed_bams", mode: "copy"

    // [sample_ID, [sample_ID.bam]]
    input:
    tuple val(sample_id), val(library), val(run), path(bam)

    output:
    path("${sample_id}-sorted-rg-md.bam")

    // Add read groups
    // Mark duplicates
    // Index bam
    script:
    """
    gatk AddOrReplaceReadGroups \
        I=${bam} \
        O=${sample_id}-sorted-rg.bam \
        RGID=${sample_id} \
        RGLB=${library} \
        RGPL=ILLUMINA \
        RGPU=${run} \
        RGSM=${sample_id}

    gatk MarkDuplicates \
        I=${sample_id}-sorted-rg.bam \
        O=${sample_id}-sorted-rg-md.bam \
        M=${sample_id}-metrics.txt

    samtools index ${sample_id}-sorted-rg-md.bam
    rm ${sample_id}-sorted-rg.bam
    """
}

process CALL_VARIANTS {

    // Directives
    cpus params.cpus
    publishDir "${params.outdir}", mode: "copy"

    input:
    path(bam)

    output:
    path("${params.vcf}.pileup")
    path("${params.vcf}.vcf")

    // Run bcftools mpileup
    // Run bcftools call
    script:
    """
    bcftools mpileup \
        --threads ${task.cpus} \
        --fasta-ref ${params.genome} \
        --min-MQ 20 \
        --min-BQ 30 \
        --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
        --output ${params.vcf}.pileup \
        --output-type z \
        ${bam}
    
    bcftools call \
        --threads ${task.cpus} \
        --multiallelic-caller \
        --variants-only \
        --ploidy 2 \
        --output-type v \
        --output ${params.vcf}.vcf \
        ${params.vcf}.pileup
    """
}
