# Trim Paired-End Reads

The [`fastp.nf`](https://github.com/Tom-Jenkins/nextflow-pipelines/blob/main/src/fastp.nf) nextflow script will take any number of samples with paired-end reads in FASTQ format and output trimmed reads using [fastp](https://github.com/OpenGene/fastp). 

## Dependencies (version tested)
* Nextflow (24.04.4)
* Java (18.0.2.1)
* Python (3.10)
* fastp (0.23.4)

## Conda Environment

Create environment using conda:   
`conda env create -f ./nextflow-pipelines/env/fastp.yml`  

Activate conda environment:  
`conda activate fastp` or `source activate fastp`

## Usage
```
#!/bin/bash

# Activate conda environment
conda activate fastp

# Variables
cpus=16
reads=/path/to/reads/directory
adapters=./nextflow-pipelines/misc/adapters.fasta
outdir=/path/to/output/directory

# Run pipeline
nextflow run ./nextflow-pipelines/src/fastp.nf \
    --reads ${reads} \
    --suffix "_{1,2}.fastq.gz" \
    --adapters ${adapters} \
    --filter "--qualified_quality_phred 30 --length_required 100 --trim_poly_g" \
    --outdir ${outdir} \
    --cpus ${cpus}
```

| Parameter&nbsp;&nbsp;&nbsp;&nbsp; | Description
| :- | :-
| `--reads` | input directory containing FASTQ files
| `--suffix` | string denoting the suffix after a sample name and read1 and read2 in the paired reads {1,2}
| `--adapters` | path to FASTA file with adapter sequences
| `--filter` | string denoting parameters passed to fastp
| `--outdir` | output directory
| `--test` | prints out a tuple of the sample ID and paths to the input paired reads (dry run)
| `--cpus` | number of cpus

**Example input:**  
```
$ ls raw_reads/
SampleID_01_1.fastq.gz SampleID_02_1.fastq.gz
SampleID_01_2.fastq.gz SampleID_02_2.fastq.gz
```

## Output

The `fastp.nf` script outputs the results to a directory called `${outdir}/trimmed_reads`. This directory contains the trimmed reads for each paired sample.

**Example output:**  
```
$ ls trimmed_reads/
SampleID_01_1.fp.fq.gz SampleID_02_1.fp.fq.gz
SampleID_01_2.fp.fq.gz SampleID_02_2.fp.fq.gz
```
