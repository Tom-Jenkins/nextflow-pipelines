#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.nano_reads = "${PWD}"
params.illumina_reads = "${PWD}"
params.pcalcareumREF = "${PWD}"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to the console
log.info """\
         A S S E M B L Y - N F   P I P E L I N E
         ===================================
         Input Nanopore directory: ${params.nano_reads}
         Input Illumina directory: ${params.illumina_reads}
         Output directory: ${params.outdir}
         Phymatolithon calcareum reference: ${params.pcalcareumREF}
         Number of threads: ${params.cpus}
         """
         .stripIndent()


// Define workflow
workflow {

    // Import Nanopore reads and store sample ID in tuple
    nanopore_reads_ch = Channel
        .fromPath("${params.nano_reads}/*.fastq.gz", checkIfExists: false).map {
            tuple ( it.name.split(".fastq.gz")[0], it)
        }
        .ifEmpty { error "No reads matching the pattern `*.fastq.gz`" }
        // .view()

    // Detect and remove adapters using Porechop v0.2.4
    porechop_ch = PORECHOP(nanopore_reads_ch)

    // Enrich for reads belonging to maerl species
    // i. Align reads to P. calcareum reference (Mor02) (Jenkins et al. 2021)
    // ii. Export mapped reads to new FASTQ file
    enrich_ch = ALIGN_TO_MAERL_REF(porechop_ch)

    // Assemble enriched maerl reads Using Flye
    flye_ch = FLYE(enrich_ch)

    // Correct assembly using Medaka (Nanopore reads)
    medaka_ch = MEDAKA(flye_ch)

    // Polish assembly using Polypolish (Illumina reads)
    polypolish_ch = POLYPOLISH(medaka_ch)

    // Map Nanopore reads (used to build assembly) to assembly
    MAP_READS(polypolish_ch)
}

// PORECHOP
process PORECHOP {

    // Directives
    cpus params.cpus
    // errorStrategy "ignore"
    conda "/lustre/home/tj311/software/mambaforge3/envs/porechop"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz")

    script:
    """
    porechop -i ${reads} -o ${sample_id}.fq -t ${params.cpus}
    gzip ${sample_id}.fq
    """
}

// ENRICHMENT
process ALIGN_TO_MAERL_REF {

    // Directives
    cpus params.cpus
    // errorStrategy "ignore"
    publishDir "${params.outdir}/${sample_id}", mode: "copy", pattern: "*.fq"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq"), optional: true

    // Align Nanopore reads to P. calcareum reference using minimap2
    // Convert to BAM and sort
    // Convert BAM to FASTQ |
    // remove unmapped reads (-F 4) |
    // remove duplicate sequences (seqkit rmdup -s) |
    // make sequence names unique (seqkit rename -n)
    script:
    """
    minimap2 -t ${params.cpus} -o enrich_${sample_id}.sam -ax map-ont ${params.pcalcareumREF} ${reads}

    samtools view -@ ${params.cpus} -b enrich_${sample_id}.sam > enrich_${sample_id}.bam

    samtools sort -@ ${params.cpus} enrich_${sample_id}.bam > enrich_${sample_id}.sort.bam
 
    rm enrich_${sample_id}.sam enrich_${sample_id}.bam
    
    samtools flagstat -@ ${params.cpus} enrich_${sample_id}.sort.bam > enrich_${sample_id}.sort.bam.stats

    samtools fastq -@ ${params.cpus} -F 4 enrich_${sample_id}.sort.bam | seqkit rmdup -s -j ${params.cpus}  | seqkit rename -n -j ${params.cpus} > enrich_${sample_id}.fq
    """
}

// FLYE ASSEMBLY
process FLYE {

    // Directives
    cpus params.cpus
    // errorStrategy "ignore"
    publishDir "${params.outdir}/${sample_id}", mode: "copy", pattern: "*"

    input:
    tuple val(sample_id), path(target_reads)

    output:
    tuple val(sample_id), path("assembly.fasta"), optional: true

    script:
    """
    flye --nano-raw ${target_reads} --out-dir . --min-overlap 5000 --threads ${params.cpus}
    """
}

// MEDAKA
process MEDAKA {

    // Directives
    cpus params.cpus
    // errorStrategy "ignore"
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(flye_assembly)

    output:
    tuple val(sample_id), path("consensus.fasta"), optional: true

    script:
    """
    medaka_consensus \
        -i ${params.outdir}/${sample_id}/enrich_${sample_id}.fq \
        -d ${flye_assembly} \
        -o . \
        -m r941_prom_sup_g507 \
        -t ${params.cpus}
    """
}

// POLYPOLISH
process POLYPOLISH {

    // Directives
    cpus params.cpus
    // errorStrategy "ignore"
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(medaka_assembly)

    output:
    tuple val(sample_id), path("polished.fasta"), optional: true

    // Index assembly
    // Align Illumina reads to genome separately (required for polypolish)
    // Execute polypolish
    script:
    """
    bwa-mem2 index ${medaka_assembly}

    bwa-mem2 mem \
        -t ${params.cpus} \
        -a ${medaka_assembly} \
        ${params.illumina_reads}/${sample_id}_trim_R1.fq.gz > alignments_1.sam
    
    bwa-mem2 mem \
        -t ${params.cpus} \
        -a ${medaka_assembly} \
        ${params.illumina_reads}/${sample_id}_trim_R2.fq.gz > alignments_2.sam

    polypolish_insert_filter.py \
        --in1 alignments_1.sam \
        --in2 alignments_2.sam \
        --out1 filtered_1.sam \
        --out2 filtered_2.sam
    
    polypolish \
        ${medaka_assembly} \
        filtered_1.sam \
        filtered_2.sam | \
        sed 's/_polypolish//' > polished.fasta
    """
}

// MAP READS
process MAP_READS {

    // Directives
    cpus params.cpus
    // errorStrategy "ignore"
    publishDir "${params.outdir}/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(polished_assembly)

    output:
    path("polished.sorted.bam")
    path("polished.sorted.bam.stats")
    path("polished.sorted.bam.*")

    // Align Nanopore reads using minimap2
    // Filter alignments using samtools
    // -F 2304: remove secondary and supplementary reads
    script:
    """
    minimap2 \
        -t ${params.cpus} \
        -o polished.sam \
        -ax map-ont ${polished_assembly} ${params.outdir}/${sample_id}/enrich_${sample_id}.fq

    samtools view -@ ${params.cpus} -F 2304 -b polished.sam > polished.bam

    samtools sort -@ ${params.cpus} polished.bam -o polished.sorted.bam

    rm polished.sam polished.bam

    samtools flagstat -@ ${params.cpus} polished.sorted.bam > polished.sorted.bam.stats

    samtools index -@ ${params.cpus} polished.sorted.bam
    """
}