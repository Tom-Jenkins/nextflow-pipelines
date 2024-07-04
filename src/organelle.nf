#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Notes:
// Make sure the seed sequences are indexed with bowtie2 prior to running pipeline
// bowtie2-build maerl-mitochondrion-seeds.fa maerl-mitochondrion-seeds.fa
// bowtie2-build maerl-chloroplast-seeds.fa maerl-chloroplast-seeds.fa

// Example usage:


// Parameters
params.reads = "${PWD}"
params.reads_suffix = "_{1,2}.fastq.gz"
params.mitogenome = false
params.plastome = false
params.mito_seed = ""
params.plastid_seed = ""
params.mito_ref = []
mito_ref = params.mito_ref ? params.mito_ref?.tokenize(",") : "None"
params.plastid_ref = []
plastid_ref = params.plastid_ref ? params.plastid_ref?.tokenize(",") : "None"
params.outdir = "${PWD}"
params.test = false
params.cpus = 16

// Print parameters to the console
log.info """\
         O R G A N E L L E - N F   P I P E L I N E
         ===================================
         Input directory: ${params.reads}
         Output directory: ${params.outdir}
         Mitochondrial seed: ${params.mito_seed != "" ? params.mito_seed : "None"}
         Mitochondrial reference(s): ${mito_ref}
         Chloroplast seed: ${params.plastid_seed != "" ? params.plastid_seed : "None"}
         Chloroplast reference(s): ${plastid_ref}
         Number of threads: ${params.cpus}
         Script version: v0.1
         """
         .stripIndent()

// Define workflow
workflow {

    // Import trimmed reads
    reads_ch = Channel
        .fromFilePairs("${params.reads}/*${params.reads_suffix}", checkIfExists: false)
        .ifEmpty { error "No paired reads found with the pattern `*${params.reads_suffix}`" }

    // Test run to view parameters and contents of reads_ch
    if ( params.test ) {
        reads_ch.view()
    }
    // Run main pipeline
    else {

        // Mitochondrial and Plastome genome assembly and annotation
        if ( params.mitogenome && params.plastome ) {
            mitoreads_ch = ALIGN_TO_MITOGENOME_SEEDS(reads_ch)
            mitogenome_ch = ASSEMBLE_MITOGENOME(mitoreads_ch)
            ANNOTATE_MITOGENOME(mitogenome_ch, mito_ref)
            plastidreads_ch = ALIGN_TO_PLASTOME_SEEDS(reads_ch)
            plastome_ch = ASSEMBLE_PLASTOME(plastidreads_ch)
            ANNOTATE_PLASTOME(plastome_ch, plastid_ref)
        }

        // Mitogenome genome assembly and annotation
        else if ( params.mitogenome ) {
            mitoreads_ch = ALIGN_TO_MITOGENOME_SEEDS(reads_ch)
            mitogenome_ch = ASSEMBLE_MITOGENOME(mitoreads_ch)
            ANNOTATE_MITOGENOME(mitogenome_ch, mito_ref)
        }
        
        // Chloroplast genome assembly and annotation
        else if ( params.plastome ) {
            plastidreads_ch = ALIGN_TO_PLASTOME_SEEDS(reads_ch)
            plastome_ch = ASSEMBLE_PLASTOME(plastidreads_ch)
            ANNOTATE_PLASTOME(plastome_ch, plastid_ref)
        }

        // Print message
        else {
            println("Error: include a --mitogenome or --plastome in the command.")
        }
    }
}


process ALIGN_TO_MITOGENOME_SEEDS {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1

    // [sample_ID, [read1.fq.gz, read2.fq.gz]]
    input:
    tuple val(sample_id), path(reads)

    // [sample_ID, sample_ID_1.fq, sample_ID_2.fq, sample_ID_unpaired.fq]]
    output:
    tuple val(sample_id), path("*.fq"), optional: true

    script:
    """
    bowtie2 -p ${task.cpus} --very-sensitive-local -x ${params.mito_seed} -1 ${reads[0]} -2 ${reads[1]} -S out.sam
    samtools view -F 260 -@ ${task.cpus} -b out.sam > out.bam
    samtools sort -n -@ ${task.cpus} out.bam > out.sorted.bam
    rm out.sam out.bam
    samtools fastq -1 ${sample_id}_1.fq -2 ${sample_id}_2.fq -s ${sample_id}_unpaired.fq out.sorted.bam
    """
}

process ASSEMBLE_MITOGENOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1
    publishDir "${params.outdir}/mitochondrial_genomes", mode: "copy"
    // conda "/lustre/home/tj311/software/miniforge3/envs/unicycler"

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
    unicycler \
        -t ${task.cpus} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -s ${reads[2]} \
        -o ${sample_id} \
        --no_rotate \
        --keep 0
    """
}

process ANNOTATE_MITOGENOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1
    publishDir "${params.outdir}/mitochondrial_genomes/${sample_id}", mode: "copy"

    // [sample_id, [sample_ID/assembly.fasta]]
    input:
    tuple val(sample_id), path(mitogenome)
    val(mito_ref)

    output:
    path("${sample_id}.cds.fasta"), optional: true

    // Construct the python command programmatically depending on the
    // number of genbank references included in the reference param
    script:
    """
    python ~/nextflow-scripts/extract_CDS.py \\
        ${mitogenome} \\
        ${mito_ref.join(' \\')} \\
        -o ${sample_id}.cds.fasta
    """
}


process ALIGN_TO_PLASTOME_SEEDS {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1

    // [sample_ID, [read1.fq.gz, read2.fq.gz]]
    input:
    tuple val(sample_id), path(reads)

    // [sample_ID, sample_ID_1.fq, sample_ID_2.fq, sample_ID_unpaired.fq]]
    output:
    tuple val(sample_id), path("*.fq"), optional: true

    script:
    """
    bowtie2 -p ${task.cpus} --very-sensitive-local -x ${params.plastid_seed} -1 ${reads[0]} -2 ${reads[1]} -S out.sam
    samtools view -F 260 -t ${task.cpus} -b out.sam > out.bam
    samtools sort -n -t ${task.cpus} out.bam > out.sorted.bam
    rm out.sam out.bam
    samtools fastq -1 ${sample_id}_1.fq -2 ${sample_id}_2.fq -s ${sample_id}_unpaired.fq out.sorted.bam
    """
}

process ASSEMBLE_PLASTOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1
    publishDir "${params.outdir}/chloroplast_genomes", mode: "copy"
    // conda "/lustre/home/tj311/software/miniforge3/envs/unicycler"

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
    unicycler \
        -t ${task.cpus} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -s ${reads[2]} \
        -o ${sample_id} \
        --no_rotate \
        --keep 0
    """
}

process ANNOTATE_PLASTOME {

    // Directives
    cpus params.cpus
    errorStrategy "ignore"
    maxForks 1 // set maximum number of parallel tasks to 1
    publishDir "${params.outdir}/chloroplast_genomes/${sample_id}", mode: "copy"

    // [sample_id, [sample_ID/assembly.fasta]]
    input:
    tuple val(sample_id), path(plastome)
    val(plastid_ref)

    output:
    path("${sample_id}.cds.fasta"), optional: true 

    // Construct the python command programmatically depending on the
    // number of genbank references included in the reference param
    script:
    """
    python ~/nextflow-scripts/extract_CDS.py \\
        ${plastome} \\
        ${plastid_ref.join(' \\')} \\
        -o ${sample_id}.cds.fasta
    """
}