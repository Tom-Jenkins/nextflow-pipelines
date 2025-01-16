import os
import gzip
from collections import defaultdict
import re
import csv
import sys

def count_reads_in_fastq(file_path):
    """
    Counts the number of reads in a FASTQ file.
    Assumes that every four lines correspond to a single read.
    """
    count = 0
    # Open a gzipped FASTQ file specified by file_path in text mode ("rt").
    # The gzip.open function decompresses the file on the fly so it can be read line by line without extracting it manually.
    # The with statement ensures that the file is properly closed after the block is executed, even if an error occurs.
    with gzip.open(file_path, "rt") as f:
    # For loop iterates over each line in the opened file.
    # The underscore _ is used as a loop variable because the actual content of the line isn't needed — only the fact that a line exists is important.
        for _ in f:
            # Count each line.
            # At the end of the loop, the count variable will represent the total number of lines in the file.
            count += 1
        # Since each read in a FASTQ file is represented by four lines, this calculates the total number of reads.
        return count // 4


def count_reads_in_directory(directory):
    """
    Counts the reads in paired FASTQ files in a directory.
    Returns a dictionary with the total reads per sample ID.
    """
    # The defaultdict function provides a default value for a key that doesn't exist in the dictionary.
    # The default value is determined by a factory function (e.g., int, list, set) provided when creating the defaultdict.
    # A regular dict requires explicit handling for missing keys, such as using 'dict.get(key, default)' or checking with 'if key in dict'.
    sample_counts = defaultdict(int)

    # Create list of files in the directory with the suffix ("fq.gz")
    files = [file for file in os.listdir(directory) if file.endswith(".fq.gz") or file.endswith(".fastq.gz") or file.endswith(".fq") or file.endswith(".fastq")]
    # return(files)

    # Regular expression to match filenames preceding and including _1 or _2
    pattern = "(.*?)_[12]"
    
    for file_name in files:
        # Match the filename pattern
        match = re.match(pattern, file_name)
        if match:
            # Extract name before _1 or _2
            sample_id = match.group(1) 

            # Get the full path of the file
            file_path = os.path.join(directory, file_name)

            # Count the read in the FASTQ file and add to the sample's total
            sample_counts[sample_id] += count_reads_in_fastq(file_path)
    
    # Write data to a CSV file
    return(sample_counts)


def export_defaultdict_as_csv(dictionary, output_csv):
    """
    Write a CSV from a default dictionary.
    """
    # Specify the output file name
    output_file = output_csv

    # Write data to a CSV file
    with open(output_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        # Header
        writer.writerow(["Sample_ID", "Number_of_reads"])
        # Data
        for key, value in dictionary.items():
            writer.writerow([key, value])


# This line checks if the script is being executed as the main program (not imported as a module in another script).
# If the script is run directly, the block of code under this if statement will execute.
# Purpose: ensures that the script’s main functionality runs only when executed directly.
if __name__ == "__main__":
    # Checks if fewer than three arguments were provided
    if len(sys.argv) < 3:
        print("Usage: python count_reads_in_fastq.py <FASTQ_directory> <output_csv>")
    else:
        input_directory = sys.argv[1]
        output_csv = sys.argv[2]
        data = count_reads_in_directory(input_directory)
        export_defaultdict_as_csv(data, output_csv)
        print(f"Completed. The output was saved in {output_csv}")
