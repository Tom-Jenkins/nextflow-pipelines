# Reference genome assembly
The [`variantcalling.nf`](https://github.com/Tom-Jenkins/maerl-wgs-pipelines/blob/main/src/variantcalling.nf) nextflow script will take any number of samples with paired-end reads in FASTQ format, map reads using Bowtie2, processed bam fildes, and finally call variants using BCFtools v1.21 and/or Freebayes v1.3.6. If part of the pipeline is unsuccessful for a sample then these errors are ignored.
Pipeline flowchart:
```mermaid
graph TD;
    FASTQ-->|bowtie2| n1[Align Reads to Reference Genome];
    n1-->|GATK4| n2[Process BAM Files];
    n2-->|| n3[Create BAM List];
    subgraph I 
        c1-->|BCFtools| c2[Call Variants];
    end
    subgraph I 
        c1-->|Freebayes| c2[Call Variants];
    end
```

## Dependencies (version tested)
* Nextflow (24.04.4)
* Java (18.0.2.1)
* Python (3.10)
* bowtie2 (2.5.3)
* SAMtools (1.19.2)
* GATK4 (4.5)
* BCFtools (1.21)
* Freebayes (1.3.6)

## Conda Environment

Create environment using conda:   
`conda env create -f ./nextflow-pipelines/env/variantcalling.yml`  

Activate conda environment:  
`conda activate variantcalling` or `source activate variantcalling`

## Setup
### Sample Sheet
The sample sheet is a __CSV__ file with five columns. Make sure the column headings are exactly the same as below.

Column 1: __sample__ (sample name)  
Column 2: __library__ (sample name if individual libraries built for each sample)  
Column 3: __run__ (if samples were sequenced on different runs use this column to specify)  
Column 4: __read1__ (absolute path to read1.fq.gz)  
Column 5: __read2__ (absolute path to read2.fq.gz)  

Example of a sample_sheet.csv:
| sample | library | run | read1 | read2
| ---  | --- | ---  | --- | ---
| sample_01 | sample_01 | Dec_2022 | /home/path/read1.fq.gz | /home/path/read2.fq.gz
| sample_02 | sample_02 | Dec_2022 | /home/path/read1.fq.gz | /home/path/read2.fq.gz
| sample_03 | sample_03 | Jun_2023 | /home/path/read1.fq.gz | /home/path/read2.fq.gz
| sample_04 | sample_04 | Jun_2023 | /home/path/read1.fq.gz | /home/path/read2.fq.gz 

## Usage
```
#!/bin/bash

# Activate conda environment
conda activate variantcalling

# Variables
cpus=20
sampleSheet=/path/to/sample/sheet.csv
genome=/path/to/reference/genome.fasta
outdir=/path/to/output/directory

# Index reference genome
bowtie2-build ${genome} ${genome}

# Run pipeline
nextflow run ~/nextflow-pipelines/src/variantcalling.nf \
    --sampleSheet ${sampleSheet} \
    --genome ${genome} \
    --outdir ${outdir} \
    --variantCaller "both" \
    --vcf "variants" \
    --cpus 20
```