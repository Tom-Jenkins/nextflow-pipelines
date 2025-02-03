#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.sampleSheet = "${PWD}/sample_sheet.csv"
params.genome = "${PWD}"
params.outdir = "${PWD}"
params.variantCaller = "both" // "bcftools" and "freebayes" or "both"
params.bcftools_mpileup = "--min-MQ 40 --min-BQ 30"
params.bcftools_call = "--ploidy 2 --multiallelic-caller --variants-only"
params.freebayes_params = "-p 2 --min-mapping-quality 40 --min-base-quality 30 --genotype-qualities"
params.vcf = "variants"
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         V A R I A N T - C A L L I N G - N F   P I P E L I N E
         ===================================
         Input reads sample sheet: ${params.sampleSheet}
         Input reference genome: ${params.genome}
         Output directory: ${params.outdir}
         Variant caller (BCFtools or Freebayes): ${params.variantCaller}
         BCFtools mpileup parameter(s): ${params.variantCaller == "bcftools" || params.variantCaller == "both" ? params.bcftools_mpileup : "n/a"}
         BCFtools call parameter(s): ${params.variantCaller == "bcftools"  || params.variantCaller == "both" ? params.bcftools_call : "n/a"}
         Freebayes parameter(s): ${params.variantCaller == "freebayes"  || params.variantCaller == "both" ? params.freebayes_params : "n/a"}
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

        // Process BAM files
        processed_bam_ch = PROCESS_BAM(bam_ch)
        
        // Collect BAM files and output text file with one BAM per line
        collected_bam_bcftools_ch = processed_bam_ch.collect()
        collected_bam_freebayes_ch = processed_bam_ch.collect()

        // Process bam files and call variants using bcftools and freebayes
        if ( params.variantCaller == "both" ) {
            CALL_VARIANTS_BCFTOOLS(collected_bam_bcftools_ch)
            CALL_VARIANTS_FREEBAYES(collected_bam_freebayes_ch)
        }

        // Process bam files and call variants using bcftools only
        if ( params.variantCaller == "bcftools" ) {
            CALL_VARIANTS_BCFTOOLS(collected_bam_bcftools_ch)
        }

        // Process bam files and call variants using freebayes only
        if ( params.variantCaller == "freebayes" ) {
            CALL_VARIANTS_FREEBAYES(collected_bam_freebayes_ch)
        }       
    }
}


process ALIGN_TO_REF_GENOME {

    // Directives
    errorStrategy "ignore"

    input:
    tuple val(sample_id), val(library), val(run), path(reads)

    // [sample_ID, [sample_ID.bam]]
    output:
    tuple val(sample_id), val(library), val(run), path("${sample_id}.bam"), optional: true

    // Align reads using bowtie2
    // Filter out unmapped reads and convert to bam
    script:
    """
    bowtie2 -p ${params.cpus} -x ${params.genome} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam
    samtools view -F 4 -@ ${params.cpus} -b ${sample_id}.sam > ${sample_id}.bam
    rm ${sample_id}.sam
    """
}

process PROCESS_BAM {

    // Directives
    errorStrategy "ignore"
    publishDir "${params.outdir}/processed_bams", mode: "copy"

    // [sample_ID, [sample_ID.bam]]
    input:
    tuple val(sample_id), val(library), val(run), path(bam)

    output:
    path("${sample_id}.rg.md.bam"), optional: true

    // Add read groups
    // Mark duplicates
    script:
    """
    echo Sorting BAM by name ...
    samtools sort -@ ${params.cpus} -n -m 4G ${bam} > ${sample_id}.sorted.temp

    echo Adding read groups ...
    samtools addreplacerg -@ ${params.cpus} -r "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tLB:${library}\\tPU:${run}\\tSM:${sample_id}" ${sample_id}.sorted.temp > ${sample_id}.rg.bam
    
    echo Adding ms and MC tags for markdup ...
    samtools fixmate -@ ${params.cpus} -m ${sample_id}.rg.bam ${sample_id}.rg.temp

    echo Sorting BAM by coordinate ...
    samtools sort -@ ${params.cpus} ${sample_id}.rg.temp > ${sample_id}.rg.sorted.temp

    echo Marking duplicates ...
    samtools markdup -@ ${params.cpus} ${sample_id}.rg.sorted.temp ${sample_id}.rg.md.bam
    
    rm ${sample_id}.sorted.temp ${sample_id}.rg.bam ${sample_id}.rg.temp ${sample_id}.rg.sorted.temp    
    """
}

process CALL_VARIANTS_BCFTOOLS {

    // Directives
    errorStrategy "ignore"
    publishDir "${params.outdir}", mode: "copy"

    input:
    path(bam)

    output:
    // path("${params.vcf}.pileup")
    path("${params.vcf}_bcftools.vcf.gz")

    // Run bcftools mpileup
    // Run bcftools call
    script:
    """
    echo ${bam} | tr " " "\n" > bamlist.txt

    bcftools mpileup \
        --threads ${params.cpus} \
        --fasta-ref ${params.genome} \
        ${params.bcftools_mpileup} \
        --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
        --output ${params.vcf}.pileup \
        --output-type z \
        --bam-list bamlist.txt
    
    bcftools call \
        --threads ${params.cpus} \
        ${params.bcftools_call} \
        --output-type v \
        --output ${params.vcf}_bcftools.vcf \
        ${params.vcf}.pileup
    
    gzip ${params.vcf}_bcftools.vcf
    """
}

process CALL_VARIANTS_FREEBAYES {

    // Directives
    errorStrategy "ignore"
    publishDir "${params.outdir}", mode: "copy"

    input:
    path(bam)

    output:
    path("${params.vcf}_freebayes.vcf.gz")

    script:
    """
    echo ${bam} | tr " " "\n" > bamlist.txt

    export TMPDIR=${params.outdir}

    samtools faidx ${params.genome}

    samtools index -M *.bam

    freebayes-parallel <(fasta_generate_regions.py ${params.genome}.fai 100000) ${params.cpus} \
        --fasta-reference ${params.genome} \
        --bam-list bamlist.txt \
        ${params.freebayes_params} \
        > ${params.vcf}_freebayes.vcf

    gzip ${params.vcf}_freebayes.vcf
    """
}
