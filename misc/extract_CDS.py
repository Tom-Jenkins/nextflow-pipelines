import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def extract_reference_genes(reference_genbanks):
    reference_genes = []

    for reference_genbank in reference_genbanks:
        genes = {}
        with open(reference_genbank, "r") as gb_fh:
            for record in SeqIO.parse(gb_fh, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        gene_seq = feature.extract(record.seq)
                        gene_id = feature.qualifiers.get("gene", ["Unknown_gene"])[0]
                        gene_name = feature.qualifiers.get("product", ["Unknown_product"])[0]
                        gene_record = SeqRecord(gene_seq, id=gene_id, description=gene_name)
                        genes[gene_id] = gene_record
        reference_genes.append(genes)

    return reference_genes

def run_blast(input_fasta, reference_fastas):
    blast_outputs = []

    for reference_fasta in reference_fastas:
        blast_output = f"blast_output_{os.path.basename(reference_fasta)}.txt"
        blast_cmd = f"blastn -query {input_fasta} -subject {reference_fasta} -out {blast_output} -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'"
        subprocess.run(blast_cmd, shell=True)
        blast_outputs.append((blast_output, reference_fasta))

    return blast_outputs

def compare_blast_outputs(blast_outputs):
    return max(blast_outputs, key=lambda x: sum(1 for line in open(x[0])))

def annotate_input_fasta(input_fasta, blast_output, output_fasta):
    annotated_sequences = {}

    output_file_name = os.path.splitext(os.path.basename(output_fasta))[0]
    
    with open(blast_output, "r") as blast_fh, open(input_fasta, "r") as input_fh:
        input_records = SeqIO.to_dict(SeqIO.parse(input_fh, "fasta"))
        for line in blast_fh:
            fields = line.strip().split("\t")
            query_id = fields[0]
            subject_id = fields[1]
            if query_id in input_records:
                query_record = input_records[query_id]
                length = int(fields[3])
                qstart = int(fields[4])
                qend = int(fields[5])
                annotated_seq = str(query_record.seq[qstart-1:qend])
                if subject_id not in annotated_sequences:
                    annotated_id = f"{subject_id} {output_file_name}"
                    annotated_sequences[subject_id] = SeqRecord(Seq(annotated_seq), id=annotated_id, description=query_id)

    with open(output_fasta, "w") as out_fh:
        SeqIO.write(annotated_sequences.values(), out_fh, "fasta")


def main(input_fasta, reference_genbanks, output_fasta):
    reference_genes = extract_reference_genes(reference_genbanks)
    reference_fastas = []
    for i, genes in enumerate(reference_genes, start=1):
        reference_fasta = f"reference_{i}.fasta"
        with open(reference_fasta, "w") as ref_fh:
            SeqIO.write(genes.values(), ref_fh, "fasta")
        reference_fastas.append(reference_fasta)

    blast_outputs = run_blast(input_fasta, reference_fastas)
    best_blast_output, best_reference_fasta = compare_blast_outputs(blast_outputs)

    annotate_input_fasta(input_fasta, best_blast_output, output_fasta)

    # Clean up intermediate files
    # for fasta in reference_fastas:
    #     os.remove(fasta)
    # for blast_output, _ in blast_outputs:
    #     os.remove(blast_output)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Annotate mitochondrial genes in input FASTA file using reference GenBank files.")
    parser.add_argument("input_fasta", help="Input FASTA file to be annotated")
    parser.add_argument("reference_genbanks", nargs="+", help="Reference GenBank files containing annotated mitochondrial genes")
    parser.add_argument("-o", "--output_fasta", default="annotated_input.fasta", help="Output annotated FASTA file (default: annotated_input.fasta)")
    args = parser.parse_args()

    main(args.input_fasta, args.reference_genbanks, args.output_fasta)
