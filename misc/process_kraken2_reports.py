import csv
import sys
import os

def process_kraken2_outputs(input_files, output_csv):
    # Prepare the data for CSV
    csv_data = []

    # Define the rank codes we are interested in
    ranks_of_interest = {'U', 'R', 'D', 'P'}

    for input_file in input_files:
        # Extract the prefix of the input file name
        prefix = os.path.splitext(os.path.basename(input_file))[0]
        
        # Open the input file and read the lines
        with open(input_file, 'r') as file:
            lines = file.readlines()
        
        # Process each line
        for line in lines:
            # Split the line into columns
            columns = line.split('\t')
            
            # Extract the relevant columns
            percentage = columns[0].strip()
            num_fragments_clade = columns[1].strip()
            rank_code = columns[3].strip()
            scientific_name = columns[5].strip()
            
            # Check if the rank code starts with one of the ranks of interest
            if rank_code[0] in ranks_of_interest:
                csv_data.append([prefix, rank_code, scientific_name, num_fragments_clade, percentage])

    # Write the data to a CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write the header
        csv_writer.writerow(['Prefix', 'Rank_code', 'Scientific_name', 'Number_of_fragments', 'Percentage'])
        # Write the data
        csv_writer.writerows(csv_data)

    print(f"Data successfully written to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <output_csv> <input_file1> <input_file2> ...")
    else:
        output_csv = sys.argv[1]
        input_files = sys.argv[2:]
        process_kraken2_outputs(input_files, output_csv)
