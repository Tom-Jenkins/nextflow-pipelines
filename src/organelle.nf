#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "${PWD}"
params.mito_seed = "${PWD}/maerl-mitochondrion-seeds.fa"
params.mito_ref = "${PWD}"
params.plastid_seed = "${PWD}/maerl-chloroplast-seeds.fa"
params.plastid_ref = "${PWD}"
params.outdir = "${PWD}"
params.cpus = 16

// Print parameters to the console
log.info """\
         O R G A N E L L E - N F   P I P E L I N E
         ===================================
         Input directory: ${params.reads}
         Output directory: ${params.outdir}
         Mitochondrial seed: ${params.mito_seed}
         Mitochondrial reference: ${params.mito_ref}
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
    ANNOTATE_MITOGENOME(mitogenome_ch)
    
}


process ALIGN_TO_GENBANK_MITOGENOMES {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"

    // [sample_ID, [sample_ID_trim_R1.fq.gz, sample_ID_trim_R2.fq.gz]]
    input:
    tuple val(sample_id), path(reads)

    // [sample_ID, sample_ID_1.fq, sample_ID_2.fq, sample_ID_unpaired.fq]]
    output:
    tuple val(sample_id), path("*.fq"), optional: true

    script:
    """
    bwa-mem2 index ${params.mito_seed}
    bwa-mem2 mem -t ${params.cpus} -o out.sam ${params.mito_seed} ${reads[0]} ${reads[1]}
    samtools view -F 4 -t ${params.cpus} -b out.sam > out.bam
    samtools sort -n -t ${params.cpus} out.bam > out.sorted.bam
    samtools fastq -1 ${sample_id}_1.fq -2 ${sample_id}_2.fq -s ${sample_id}_unpaired.fq out.sorted.bam
    """
}

process ASSEMBLE_MITOGENOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    publishDir "${params.outdir}/mitochondrial_genomes", mode: "copy"
    conda "/lustre/home/tj311/software/mambaforge3/envs/unicycler"

    // [sample_ID, sample_ID_1.fq, sample_ID_2.fq, sample_ID_unpaired.fq]]
    input:
    tuple val(sample_id), path(reads)

    // [sample_id, [sample_ID/assembly.fasta]]
    output:
    tuple val(sample_id), path("${sample_id}/assembly.fasta"), optional: true
    // tuple val(sample_id), path("${sample_id}/assembly.gfa"), optional: true
    // tuple val(sample_id), path("${sample_id}/assembly.log"), optional: true

    script:
    """
    unicycler -t ${params.cpus} -1 ${reads[0]} -2 ${reads[1]} -s ${reads[2]} -o ${sample_id} --keep 0
    """
}

process ANNOTATE_MITOGENOME {

    // Directives
    cpus params.cpus
    publishDir "${params.outdir}/mitochondrial_genomes/${sample_id}/annotation", mode: "copy"

    // [sample_id, [sample_ID/assembly.fasta]]
    input:
    tuple val(sample_id), path(mitogenome)

    // [sample_id_MitoFinder.log]
    // [sample_id/*_Results/*all_files]
    output:
    path("*.log"), optional: true 
    path("${sample_id}/*_Results/*"), optional: true

    script:
    """
    mitofinder \
        --seqid ${sample_id} \
        --assembly ${mitogenome} \
        --refseq ${params.mito_ref} \
        --organism 4 \
        --processors ${task.cpus} \
        --max-contig-size 30000 \
        --new-genes \
        --adjust-direction
    """
}



