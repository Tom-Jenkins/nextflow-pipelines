import os
import sys
import csv
import argparse


def paths_to_quant_files(path, salmon=False, kallisto=False):
    """
    Parse paths to Salmon or Kallisto quantification files
    Returns a dictionary in the format {sample_name: sample_path}
    """
    if salmon:
        # List of samples for Salmon output
        salmon_samples = os.listdir(path)
        
        # Paths to quant.sf files
        salmon_paths = [f"{path}/{file}/quant.sf" for file in salmon_samples]
        
        # Dictionary
        salmon_dict = {}
        for i, sample in enumerate(salmon_samples):
            salmon_dict[sample] = salmon_paths[i]
        # print(salmon_dict)
        return(salmon_dict)

    if kallisto:
        # List of samples for Kallisto output
        kallisto_samples = os.listdir(path)
        
        # Paths to abundance.tsv files
        kallisto_paths = [f"{path}/{file}/abundance.tsv" for file in kallisto_samples]

        # Dictionary
        kallisto_dict = {}
        for i, sample in enumerate(kallisto_samples):
            kallisto_dict[sample] = kallisto_paths[i]
        # print(kallisto_dict)
        return(kallisto_dict)


def extract_gene_id(dict, col=0):
    """
    Extract the gene ID from a text file, typically from column 1 in Salmon and Kallisto
    Returns a list of strings
    """
    # Extract the path for the first file from the dictionary
    file = list(dict.values())[0]

    # Open and read the file
    with open(file, "r") as f:
        lines = f.readlines()

    # Subset column 1
    geneID = []
    for line in lines[1:]: # skip the header
        columns = line.strip().split() # split by whitespace
        geneID.append(columns[col])
    # print(geneID[:3])
    return(geneID)


def merge_quant_files(dict, geneIDs, program):
    """
    Merge count data for all samples into a single file
    Returns a dictionary where the first key is the geneID and the remaining keys are the samples
    """
    # Sample names
    sample_names = list(dict.keys())
    # print(sample_names)

    # Sample file paths
    sample_paths = list(dict.values())
    # print(sample_paths)

    # Create a dictionary to store count data
    data_dict = {
        "GeneID": geneIDs,
    }
    # print(data_dict["Gene"][:3])

    # In Salmon quant.sf column 5 are counts
    if program == "salmon":
        count_column = 4

    # In Kallisto abundance.tsv column 4 are counts
    if program == "kallisto":
        count_column = 3

    # Open and read each file in sample_paths
    for i, file in enumerate(sample_paths):
        with open(file, "r") as f:
            lines = f.readlines()

            # Subset column 5
            counts = []
            for line in lines[1:]: # skip the header
                columns = line.strip().split() # split by whitespace
                counts.append(columns[count_column]) # extract counts from column
            
            # Append counts to dictionary
            data_dict[sample_names[i]] = counts
    
    # Return dictionary
    # print(data_dict.keys())
    return(data_dict)


def export_dict_as_csv(data_dict, output_csv):
    """
    Write a CSV from a dictionary
    """
    # Output file name
    output_file = output_csv

    # Transpose dictionary into rows
    rows = zip(*data.values())

    # Write data to CSV file
    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        # Header
        writer.writerow(data_dict.keys())
        # Data
        writer.writerows(rows)


# This line checks if the script is being executed as the main program (not imported as a module in another script).
# If the script is run directly, the block of code under this if statement will execute.
# Purpose: ensures that the script’s main functionality runs only when executed directly.
if __name__ == "__main__":
 
    # Parse arguments on command line
    parser = argparse.ArgumentParser(description="Merge quantification files from Salmon or Kallisto.")

    # Ensure user can only specify --salmon or --kallisto
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--salmon", action="store_true", help="Use Salmon results files")
    group.add_argument("--kallisto", action="store_true", help="Use Kallisto results files")

    # Other arguments
    parser.add_argument("results_dir", help="Path to results directory")
    parser.add_argument("output_csv", help="Path to output CSV file")

    args = parser.parse_args()

    # Determine the program chosen on command line
    if args.salmon:
        program = "salmon"
    elif args.kallisto:
        program = "kallisto"
    
    # Execute Salmon
    if program == "salmon":
        quant_dict = paths_to_quant_files(args.results_dir, salmon=True, kallisto=False)
        geneIDs = extract_gene_id(quant_dict)
        data = merge_quant_files(quant_dict, geneIDs, program=program)
        export_dict_as_csv(data, args.output_csv)
        print(f"Completed. The output was saved in {args.output_csv}")

    # Execute Kallisto
    if program == "kallisto":
        quant_dict = paths_to_quant_files(args.results_dir, salmon=False, kallisto=True)
        geneIDs = extract_gene_id(quant_dict)
        data = merge_quant_files(quant_dict, geneIDs, program=program)
        export_dict_as_csv(data, args.output_csv)
        print(f"Completed. The output was saved in {args.output_csv}")
