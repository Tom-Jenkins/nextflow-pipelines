# Quantify Expression Using Salmon or Kallisto

The [`quantifyexpression.nf`](https://github.com/Tom-Jenkins/nextflow-pipelines/blob/main/src/quantifyexpression.nf) nextflow script will take any number of paired reads in FASTQ format and quantify expression using Salmon or/and Kallisto. 

## Dependencies (version tested)
* Nextflow (24.04.4)
* Java (18.0.2.1)
* Python (3.10)
* Salmon (1.10.3)
* Kallisto (0.51.1)

## Conda Environment

Create environment using conda:   
`conda env create -f ./nextflow-pipelines/env/quantifyexpression.yml`  

Activate conda environment:  
`conda activate quantifyexpression` or `source activate quantifyexpression`

## Usage
```
#!/bin/bash

# Activate conda environment
conda activate quantifyexpression

# Variables
cpus=20
reads=/path/to/reads/directory
transcriptome=/path/to/transcriptome.fasta
outdir=/path/to/output/directory

# Run pipeline
nextflow run ~/nextflow-pipelines/src/quantifyexpression.nf \
    --reads ${reads} \
    --suffix "_{1,2}.fq.gz" \
    --transcriptome ${transcriptome} \
    --salmon \
    --salmonParams "--libType A --gcBias --seqBias" \
    --kallisto \
    --outdir ${outdir} \
    --cpus ${cpus}

# Merge quantification files
python ~/nextflow-pipelines/misc/merge_quant_files.py --salmon salmon_output/ salmon_quant_results.csv
python ~/nextflow-pipelines/misc/merge_quant_files.py --kallisto kallisto_output/ kallisto_quant_results.csv
```

| Parameters&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Description
| :--- | :---
| `--reads` | path to input directory containing FASTQ files
| `--suffix` | string denoting the suffix after a sample name and the forward (read1) and reverse (read2) designation (e.g. for read pair `sampleID_1.fq.gz` and `sampleID_2.fq.gz` set the parameter to `--suffix "_{1,2}.fq.gz"`. The name of this output directory will be called `sampleID`) 
| `--transcriptome` | path to transcriptome in FASTA format
| `--salmon` | run the pipeline using Salmon
| `--salmonParams` | string of additional arguments passed to `salmon quant` command (default: "--libType A --gcBias --seqBias")
| `--kallisto` | run the pipeline using Kallisto
| `--kallistoParams` | string of additional arguments passed to `kallisto quant` (default: " ")
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
SampleID_01/ SampleID_02/
```
