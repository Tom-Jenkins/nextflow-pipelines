# Nextflow Pipelines for Bioinformatics Analysis of Sequencing Data from Multiple Samples

> Note: this repository will be continually updated as scripts are tested. Pipelines not ready for deployment are noted as *in development* below. 

The main purpose of this repository is to provide nextflow scripts that process sequencing data from multiple samples. These scripts have been designed for a specific purpose, such as trimming reads from multiple paired samples, and are an alternative to setting up bash loops or writing other custom scripts.

## Setup
There are two ways to install the software needed to run the nextflow scripts. The recommended way is to set up a conda environment which will install all the software needed, including nextflow, to run a pipeline. Instructions for setting up conda environments are detailed in the documentation. An alternative is to install all the software needed manually and then call a nextflow script. If using this manual approach, the software dependencies for each pipeline is detailed in the documentation.

### Conda or Mamba
Install conda from [miniforge3](https://github.com/conda-forge/miniforge?tab=readme-ov-file#miniforge3) or from [miniconda3](https://docs.anaconda.com/free/miniconda).

Add conda channels:  
`conda config --add channels conda-forge`  
`conda config --add channels bioconda`

### Nextflow
If you are manually installing software, nextflow can be installed by following the installation [docs](https://www.nextflow.io/docs/latest/getstarted.html). The scripts have currently been tested with nextflow v23.10.1 and java v18.0.2.1 2022-08-18.

## Documentation

### Nextflow Pipeline (*in development*)
- [Trim Illumina paired-end reads](https://github.com/Tom-Jenkins/maerl-wgs-pipelines/blob/main/docs/01-trim-illumina-reads.md)
- [Organelle assembly and annotation](https://github.com/Tom-Jenkins/maerl-wgs-pipelines/blob/main/docs/02-organelle-assembly-annotation.md)
- [Reference genome assembly](https://github.com/Tom-Jenkins/maerl-wgs-pipelines/blob/main/docs/03-reference-genome-assembly.md)


