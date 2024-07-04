# Trim Illumina Paired-End Reads

The [`fastp.nf`](https://github.com/Tom-Jenkins/nextflow-pipelines/blob/main/src/fastp.nf) nextflow script will take any number of samples with paired-end reads in FASTQ format and output trimmed reads using [fastp](https://github.com/OpenGene/fastp). 

## Conda Environment

Create environment using conda:   
`conda env create -f ./nextflow-pipelines/env/fastp.yml`  

Create environment using mamba (faster):  
`mamba env create -f ./nextflow-pipelines/env/fastp.yml`

Activate conda environment:  
`mamba activate fastp` or `conda activate fastp` or `source activate fastp`

## Usage
```
#!/bin/bash

# Activate conda environment
mamba activate fastp

# Variables
cpus=16
reads=/path/to/reads/directory
esf=/path/to/reads/directory
adapters=./nextflow-pipelines/misc/adapters.fasta
outdir=/path/to/output/directory

# Run pipeline
nextflow run ./nextflow-pipelines/src/fastp.nf \
    --reads ${reads} \
    --reads_suffix "_{1,2}.fastq.gz" \
    --esf ${esf} \
    --esf_prefix "11171_|11002_|10628_" \
    --esf_suffix "_R{1,2}_001.fastq.gz" \
    --adapters ${adapters} \
    --outdir ${outdir} \
    --cpus ${cpus}
```

| Parameter | Description
| :- | :-
| `--reads` | input directory containing FASTQ files
| `--reads_suffix` | string denoting the suffix after a sample name and read1 and read2 in the paired reads {1,2}
| `--esf` | input directory containing FASTQ files with prefixes
| `--esf_prefix` | string denoting the prefix before a sample name
| `--esf_suffix` | string denoting the suffix after a sample name and read1 and read2 in the paired reads {1,2}
| `--adapters` | path to FASTA file with adapter sequences
| `--outdir` | output directory
| `--test` | prints out a tuple of the sample ID and paths to the input paired reads (dry run)
| `--cpus` | number of cpus
> `--reads` or `--esf` is required. Strings for `--esf_prefix` can contain a pipe `|` when multiple prefixes are present.

The `--esf` and related parameters work for sequencing files in the following format: `10628_Sample_ID_S1_R1_001.fastq.gz`.

**Example input:**  
```
$ ls raw_reads/
SampleID_01_1.fastq.gz SampleID_02_1.fastq.gz
SampleID_01_2.fastq.gz SampleID_02_2.fastq.gz
```
```
$ ls raw_reads_esf/
10628_Sample_ID_S3_R1_001.fastq.gz 10628_Sample_ID_S40_R1_001.fastq.gz
10628_Sample_ID_S3_R2_001.fastq.gz 10628_Sample_ID_S40_R2_001.fastq.gz
```

## Output

The `fastp.nf` script outputs the results to a directory called `${outdir}/trimmed_reads`. This directory contains the trimmed reads for each paired sample.

**Example output:**  
```
$ ls trimmed_reads/
SampleID_01_1.fp.fq.gz SampleID_01_1.fp.fq.gz
SampleID_01_2.fp.fq.gz SampleID_01_2.fp.fq.gz
```
