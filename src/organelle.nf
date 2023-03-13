#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.mito_seed = "${PWD}/maerl-mitochondrion-seeds.fa"
params.plastid_seed = "${PWD}/maerl-chloroplast-seeds.fa"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to the console
log.info """\
         O R G A N E L L E - N F   P I P E L I N E
         ===================================
         Input directory: ${params.reads}
         Output directory: ${params.outdir}/organelle_genomes
         Mitochondrial seeds directory: ${params.mito_seed}
         Number of threads: ${params.cpus}
         """
         .stripIndent()


// Define workflow
workflow {

    // Import reads
    reads_ch = Channel
        .fromFilePairs("${params.reads}/*_trim_R{1,2}.fq.gz", checkIfExists: false)
        .ifEmpty { error "No paired reads matching the pattern `*_trim_R{1,2}.fq.gz`" }
        // .view()

    // Mitochondrial genome assembly and annotation
    mitoreads_ch = ALIGN_TO_GENBANK_MITOGENOMES(reads_ch)
    mitogenome_ch = ASSEMBLE_MITOGENOME(mitoreads_ch)
    
}


process ALIGN_TO_GENBANK_MITOGENOMES {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_{1,2}.fq") optional true

    script:
    """
    bwa-mem2 index ${params.mito_seed}
    bwa-mem2 mem -t ${params.cpus} -o out.sam ${params.mito_seed} ${reads[0]} ${reads[1]}
    samtools view -f 2 -t ${params.cpus} -b out.sam > out.bam
    samtools sort -t ${params.cpus} out.bam > out.sorted.bam
    samtools fastq -1 ${sample_id}_1.fq -2 ${sample_id}_2.fq out.sorted.bam
    """
}

process ASSEMBLE_MITOGENOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    publishDir "${params.outdir}/mitogenome", mode: "copy"
    conda "/lustre/home/tj311/software/mambaforge3/envs/unicycler"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.fasta") optional true
    // path("${sample_id}/*.gfa") optional true
    // path("${sample_id}/*.log") optional true

    script:
    """
    unicycler -t ${params.cpus} -1 ${reads[0]} -2 ${reads[1]} -o . --keep 0
    seqkit seq --min-len 25000 --max-len 30000 assembly.fasta > ${sample_id}.fasta
    """
}

// process ANNOTATE_MITOGENOME {

//     // Directives
//     cpus params.cpus
//     publishDir "${params.outdir}/mitogenome", mode: "copy"

//     input:
//     path(mitogenome)

//     output:
//     path("${sample_id}.XXX") optional true

//     script:
//     """
//     """
// }




// // EXTRACT MITOGENOME FROM READS
// process GET_MITOGENOME {

//     // Directives
//     cpus params.cpus
//     publishDir "${params.outdir}/organelle_genomes/mitogenome", mode: "copy" // pattern: "${sample_id}/*.paired.fq"
//     conda "/lustre/home/tj311/software/mambaforge3/envs/getorganelle"

//     input:
//     tuple val(sample_id), path(reads)

//     output:
//     path("${sample_id}/*.gfa") optional true
//     path("${sample_id}/*.fasta") optional true
//     path("${sample_id}/*.txt") optional true

//     script:
//     """
//     get_organelle_from_reads.py \
//         -1 ${reads[0]} \
//         -2 ${reads[1]} \
//         -o ${sample_id} \
//         -F embplant_mt \
//         -s ${params.mito_seed} \
//         --target-genome-size 30000 \
//         --expected-min-size 25000 \
//         --expected-max-size 30000 \
//         --disentangle-time-limit 600 \
//         --max-paths-num 1 \
//         --random-seed 123 \
//         -t ${task.cpus}
//     """
// }