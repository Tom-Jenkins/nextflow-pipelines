# Trim Illumina paired-end reads

The `fastp.nf` nextflow script will take any number of samples with paired-end reads in FASTQ format and output trimmed reads using [fastp](https://github.com/OpenGene/fastp). The names of the files sent by the Sequencing Facility were in the following format: `10628_SampleID_PlateID_R1_001.fastq.gz`. Therefore, the `fastp.nf` script recognises this pattern (`*_R{1,2}_001.fastq.gz`) for each sample pair and during processing removes the project ID `10628_` from the beginning and the plate ID (e.g. `S72`) from the middle of each sample name. If you want to adapt this script to suit your own file formats, edit both the pattern recognition and the project/plate ID code segments.

## 1. Install fastp

First, create a conda environment for fastp v0.23.2.
```
# Create conda env
mamba create -n fastp -c bioconda fastp=0.23.2

# Print env paths on system
mamba env list
```
Second, edit path to the fastp conda environment in the `fastp.nf` script (line 22).
```
conda "/path/to/mambaforge3/envs/fastp"
```

## 2. Run fastp nextflow script

```
nextflow run ./src/fastp.nf --cpus 16 --reads /path/to/raw_reads/ --outdir /path/to/output_dir/
```
| Parameter | Description
| :- | :-
| --cpus | number of threads
| --reads | directory path containing input FASTQ files
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

