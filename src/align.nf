#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.suffix = "_{1,2}.fastq.gz"
params.genome = "${PWD}"
params.outdir = "${PWD}"
params.aligner = "bowtie2"
params.bowtie2 = ""
params.bwamem2 = ""
params.filter = "-F 4"
params.indexbam = false
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         A L I G N - N F   S C R I P T
         ===================================
         Input reads directory: ${params.reads}
         Input reads suffix: ${params.suffix}
         Input reference genome: ${params.genome}
         Aligner software chosen: ${params.aligner}
         Aligner parameter(s): ${params.aligner == "bowtie2" ? params.bowtie2 : params.bwamem2}
         Filtering parameter(s): ${params.filter}
         Output directory: ${params.outdir}
         Number of threads: ${params.cpus}
         Script version: v0.3
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

        // BOWTIE2
        if (params.aligner == "bowtie2") {

            // Compute and process alignments
            alignments = BOWTIE2(reads_ch)
            processed_bams = SAMTOOLS(alignments)

            // Index BAMs if parameter is set
            if (params.indexbam) INDEX_BAM(processed_bams)
        }

        // BWA-MEM2
        if (params.aligner == "bwamem2") {

            // Compute and process alignments
            alignments = BWAMEM2(reads_ch)
            processed_bams = SAMTOOLS(alignments)

            // Index BAMs if parameter is set
            if (params.indexbam) INDEX_BAM(processed_bams)
        }
    }
}

// BOWTIE2
process BOWTIE2 {

    // Directives
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1

    // [sample_ID, [read1, read2]]
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bowtie2 -p ${params.cpus} ${params.bowtie2} -x ${params.genome} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam
    samtools view -@ ${params.cpus} -b ${sample_id}.sam > ${sample_id}.bam
    rm ${sample_id}.sam
    """
}

// BWA-MEM2
process BWAMEM2 {
    
    // Directives
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1

    // [sample_ID, [read1, read2]]
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa-mem2 mem -t ${params.cpus} ${params.bwamem2} -a ${params.genome} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    samtools view -@ ${params.cpus} -b ${sample_id}.sam > ${sample_id}.bam
    rm ${sample_id}.sam
    """
}

// SAMTOOLS
process SAMTOOLS {
    
    // Directives
    publishDir "${params.outdir}", mode: "copy"
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1

    input:
    tuple val(sample_id), path(bam)

    output:
    path("${sample_id}.bam")
    // path("${sample_id}.bam.stats")

    // Default is to filter out unmapped reads (-F 4)
    script:
    """
    samtools view -@ ${params.cpus} ${params.filter} -b ${bam} > ${sample_id}.unsorted.bam
    samtools sort -@ ${params.cpus} ${sample_id}.unsorted.bam -o ${sample_id}.bam
    rm ${sample_id}.unsorted.bam
    """
}

// INDEX BAM
process INDEX_BAM {

    // Directives
    publishDir "${params.outdir}", mode: "copy"
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1

    input:
    path(bam)

    output:
    path("*.bai")

    script:
    """
    samtools index -@ ${params.cpus} --bai ${bam}
    """
}
