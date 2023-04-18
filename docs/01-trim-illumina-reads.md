# Trim Illumina paired-end reads

The [`fastp.nf`](https://github.com/Tom-Jenkins/maerl-wgs-pipelines/blob/main/src/fastp.nf) nextflow script will take any number of samples with paired-end reads in FASTQ format and output trimmed reads using [fastp](https://github.com/OpenGene/fastp). The names of the files sent by the Sequencing Facility were in the following format: `10628_Sample_ID_PlateID_R1_001.fastq.gz`. Therefore, the `fastp.nf` script recognises this pattern (`*_R{1,2}_001.fastq.gz`) for each sample pair and during processing removes the project ID `10628_` from the beginning and the plate ID (e.g. `S72`) from the middle of each sample name. This script also optionally accepts reads downloaded from the Sequence Read Archive (SRA). If you want to adapt this script to suit your own file formats, edit both the pattern recognition and the project/plate ID code segments.

**Download SRA reads for Project [PRJNA682082](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA682082)**
```
mamba create -n fastq-dl -c conda-forge -c bioconda fastq-dl=2.0.0
mamba activate fastq-dl
fastq-dl --outdir sra_reads/ --cpus 4 --provider SRA --accession PRJNA682082
```

**Example input:**
```
$ ls raw_reads
```
```
10628_Sample_ID_S3_R1_001.fastq.gz 10628_Sample_ID_S40_R1_001.fastq.gz
10628_Sample_ID_S3_R2_001.fastq.gz 10628_Sample_ID_S40_R2_001.fastq.gz
```
```
$ ls sra_reads
```
```
SRR13356831_1.fastq.gz SRR13356832_1.fastq.gz
SRR13356831_2.fastq.gz SRR13356832_2.fastq.gz
```

## 1. Install fastp

First, create a conda environment for fastp v0.23.2.
```
# Create conda env
mamba create -n fastp -c bioconda fastp=0.23.2

# Print env paths on system
mamba env list
```
Second, edit path to the fastp conda environment in the `fastp.nf` script (line 41).
```
conda "/path/to/mambaforge3/envs/fastp"
```

## 2. Run fastp nextflow script

```
nextflow run ./src/fastp.nf --cpus 16 --reads /path/to/raw_reads/ --sra /path/to/sra_reads/ --outdir /path/to/output_dir/
```
| Parameter | Description
| :- | :-
| --cpus | number of threads
| --reads | directory path containing input FASTQ files
| --sra | directory path containing input SRA FASTQ files (optional)
| --outdir | directory path for trimmed output FASTQ files

## Output

The output of `fastp.nf` is a directory called `trimmed_reads/` that is automatically created in the `--outdir` path. In this directory are the trimmed reads for each sample pair that were present in the input `--reads` directory.

**Example output:**
```
$ ls trimmed_reads
```
```
SampleID1_trim.R1.fq.gz SampleID2_trim.R1.fq.gz
SampleID1_trim.R2.fq.gz SampleID2_trim.R2.fq.gz
```

