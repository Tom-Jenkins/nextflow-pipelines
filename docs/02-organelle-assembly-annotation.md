# Organelle assembly and annotation

Conda:

- Unicycler

Additional Dependencies:

- bwa-mem2
- samtools
- MitoFinder
- PGA
- gbseqextractor
- SeqKit

**Download mitochondrial genomes from GenBank to use as seed sequences**

`efetch -db nucleotide -id MW357900.2,MW357901.2,MH281621.1 -format fasta > maerl-mitochondrion-seeds.fa`

**Download annotation files from GenBank to use as references**

`efetch -db nucleotide -id MW357900.2 -format gb > MW357900-mitochondrion-reference.gb`


**Download chloroplast genomes from GenBank to use as seed sequences**

`efetch -db nucleotide -id OQ417768.1,OQ417769.1,MH281627.1 -format fasta > maerl-chloroplast-seeds.fa`

**Download annotation files from GenBank to use as references**

`efetch -db nucleotide -id OQ417768.1 -format gb > chloroplast_references/OQ417768-chloroplast-reference.gb`
`efetch -db nucleotide -id OQ417769.1 -format gb > chloroplast_references/OQ417769-chloroplast-reference.gb`
`efetch -db nucleotide -id MH281627.1 -format gb > chloroplast_references/MH281627-chloroplast-reference.gb`

https://pypi.org/project/gbseqextractor/

cat ~/maerl/novaseq/organelle_genomes/mitochondrial_genomes/*/*/*/*Results/*final_genes_NT.fasta | seqkit grep -r -p "COX1"

cat ~/maerl/novaseq/organelle_genomes/chloroplast_genomes/*/annotation/*.cds.fasta | seqkit grep -r -p "psbA"
