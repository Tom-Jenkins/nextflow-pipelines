# Align Paired-End Reads To Reference

The [`align.nf`](https://github.com/Tom-Jenkins/nextflow-pipelines/blob/main/src/align.nf) nextflow script will take any number of paired reads in FASTQ format and output sorted BAM alignments reads using [bowtie2](https://github.com/BenLangmead/bowtie2) or [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2). 

## Dependencies (version tested)
* Nextflow (24.04.4)
* Java (18.0.2.1)
* Python (3.10)
* bowtie2 (2.5.3)
* bwa-mem2 (2.2.1)
* SAMtools (1.19.2)

## Conda Environment

Create environment using conda:   
`conda env create -f ./nextflow-pipelines/env/align.yml`  

Activate conda environment:  
`conda activate align` or `source activate align`

## Usage
```
#!/bin/bash

# Activate conda environment
conda activate align

# Variables
cpus=20
reads=/path/to/reads/directory
genome=/path/to/reference/filename.fasta
outdir=/path/to/output/directory

# Index genome
# bowtie2-build ${genome} ${genome}
# bwa-mem2 index ${genome}

# Run pipeline
nextflow run ~/nextflow-pipelines/src/align.nf \
    --reads ${reads} \
    --suffix "_{1,2}.fq.gz" \
    --reference ${genome} \
    --outdir ${outdir} \
    --aligner "bowtie2" \
    --filter "-F 4" \
    --cpus ${cpus}
```

| Parameters&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description
| :--- | :---
| `--reads` | path to input directory containing FASTQ files
| `--suffix` | string denoting the suffix after a sample name and the forward (read1) and reverse (read2) designation (e.g. for read pair `sample_1.fq.gz` and `sample_2.fq.gz` set the parameter to `--suffix "_{1,2}.fq.gz"`. The name of this BAM file will be called `sample.bam`) 
| `--reference` | path to reference FASTA file (e.g. reference genome)
| `--outdir` | path to output directory
| `--aligner` | string denoting whether to use bowtie2 `"bowtie2"` or bwa-mem2 `"bwa-mem2"` for alignments (default: `"bowtie2"`)
| `--bowtie2`| string of additional arguments passed to bowtie2 (e.g. `--bowtie2 "--sensitive --seed 123"`)
| `--bwamem2`| string of additional arguments passed to bwa-mem2 (e.g. `--bwamem2 "-M -k 19"`)
| `--filter` | string of additional arguments passed to samtools to filter BAMs (e.g. to keep only primary alignments: `--filter "-F 260"`, default is to remove unmapped reads: `--filter "-F 4"`)
| `--indexbam` | index BAMs using `samtools index` (optional)
| `--test` | prints out a tuple of the sample ID and paths to the input paired reads (dry run)
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
SampleID_01.bam SampleID_02.bam
```
