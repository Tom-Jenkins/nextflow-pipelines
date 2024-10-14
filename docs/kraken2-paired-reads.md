# Classify Paired-End Reads Using Kraken2

The [`kraken2.nf`](https://github.com/Tom-Jenkins/nextflow-pipelines/blob/main/src/kraken2.nf) nextflow script will take any number of paired reads in FASTQ format and classify reads using Kraken2. 

## Dependencies (version tested)
* Nextflow (23.10.1)
* Java (18.0.2.1)
* Python (3.10)
* Kraken2 (2.1.3)

## Conda Environment

Create environment using conda:   
`conda env create -f ./nextflow-pipelines/env/kraken2.yml`  

Create environment using mamba (faster):  
`mamba env create -f ./nextflow-pipelines/env/kraken2.yml`

Activate conda environment:  
`mamba activate kraken2` or `conda activate kraken2` or `source activate kraken2`

## Usage
```
#!/bin/bash

# Activate conda environment
mamba activate kraken2

# Variables
cpus=20
reads=/path/to/reads/directory
krakenDB=/path/to/krakenDB/
outdir=/path/to/output/directory

# Run pipeline
nextflow run ~/nextflow-pipelines/src/kraken2.nf \
    --reads ${reads} \
    --suffix "_{1,2}.fq.gz" \
    --krakenDB ${krakenDB} \
    --outdir ${outdir} \
    --cpus ${cpus}
```

| Parameters&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description
| :--- | :---
| `--reads` | path to input directory containing FASTQ files
| `--suffix` | string denoting the suffix after a sample name and the forward (read1) and reverse (read2) designation (e.g. for read pair `sample_1.fq.gz` and `sample_2.fq.gz` set the parameter to `--suffix "_{1,2}.fq.gz"`. The name of this BAM file will be called `sample.bam`) 
| `--krakenDB` | path to Kraken2 database
| `--kraken2` | string of additional arguments passed to bowtie2 (e.g. "--gzip-compressed --minimum-hit-groups 3")
| `--outdir` | path to output directory
| `--cpus` | integer denoting the number of cpus (default: `16`)


## Input

```
$ ls input/
SampleID_01_1.fq.gz SampleID_02_1.fq.gz
SampleID_01_2.fq.gz SampleID_02_2.fq.gz
```

## Output

```
$ ls output/
SampleID_01.kraken SampleID_02.kraken
SampleID_01.report SampleID_02.report
```
