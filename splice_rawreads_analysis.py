import pysam
from sys import argv
from tabulate import tabulate
import csv
import os
import pandas as pd
from tqdm import tqdm
from colorama import Fore, Style
import argparse

def calculate_splicing_efficiency(bam_file):
    # Initialize variables to keep track of reads and coverage
    total_reads = 0
    sample_name = os.path.basename(bam_file).split("_")[1].split("-")[0]

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        all_reads = list(bam.fetch())

        # find introns in all reads
        all_introns = bam.find_introns(all_reads)

        # extract only the keys, not include counts
        intron_regions = [(k[0], k[1]) for k in all_introns.keys()]

        # Initialize intron dictionary and other counters
        myintrons = {k: {"splice_support": 0, "fiveprime_unspliced": 0,"threeprime_unspliced":0,"fiveprime_boundary":0,
                         "threeprime_boundary":0} for k in intron_regions}
        splice_support = 0
        non_splice_support = 0

        # Loop through each read in the BAM file
        # Bin reads into exons or introns depending on they cigar string attributes
        for read in all_reads:
            # Initialize variablesa
            exonic_length = 0
            spliced_length = 0

            total_reads += 1

            for op, length in read.cigartuples:
                if op == 0 or op == 7 or op == 8:
                    exonic_length += length
                elif op == 3:
                    spliced_length += length

            # Calculates the reads that are supporting splice
            # Uses the splice length variable which contains the cigar strings with "N"
            if spliced_length > 0:
                running_over = [t for t in intron_regions if
                                read.reference_start < t[0] + 3 and read.reference_end > t[1] - 3]
                if len(running_over) > 0:
                    for t in running_over:
                        intron_length = t[1] - t[0]
                        if spliced_length <= intron_length <= spliced_length:
                            myintrons[t]["splice_support"] += 1
                    splice_support += 1

            # Calculates the five and threeprime unspliced reads and boundary reads
            # Boundary reads are taken into account as it's not certain if they're spliced/unspliced
            if spliced_length == 0:
                fiveprime_unspliced = [t for t in intron_regions if
                                       read.reference_start < t[0] and
                                       read.reference_end > t[0] and
                                       read.reference_end < t[1] and
                                       "N" not in read.cigarstring]

                threeprime_unspliced = [t for t in intron_regions if
                                        read.reference_start < t[1] and
                                        read.reference_end > t[1] and
                                        "N" not in read.cigarstring]

                if len(fiveprime_unspliced) > 0:
                    for t in fiveprime_unspliced:
                        myintrons[t]["fiveprime_unspliced"] += 1

                if len(threeprime_unspliced) > 0:
                    for t in threeprime_unspliced:
                        myintrons[t]["threeprime_unspliced"] += 1

                five_boundary = [t for t in intron_regions if
                                read.reference_start <= t[0] and
                                t[0] + 1 <= read.reference_end <= t[0] + 5 and
                                "N" not in read.cigarstring]
                if len(five_boundary) > 0:
                    for t in five_boundary:
                        myintrons[t]["fiveprime_boundary"] += 1

                three_boundary = [t for t in intron_regions if
                                t[1] - 5 <= read.reference_start <= t[1] - 1 and
                                read.reference_end >= t[1] and
                                "N" not in read.cigarstring]

                if len(three_boundary) > 0:
                    for t in three_boundary:
                        myintrons[t]["threeprime_boundary"] += 1

        # Print the myintrons dictionary for this BAM file
        print(f"Myintrons for {sample_name}:")
        print(myintrons)

    return myintrons
# Define the input directory for BAM files
bam_files_directory = "/homes/biertank/ella/ref_genomes/all_genomes/alignment_output_3/bamfiles/hpv16_bam"

# Get a list of BAM files in the specified directory
bam_files = [f for f in os.listdir(bam_files_directory) if f.endswith(".bam")]

# Initialize a dictionary to store raw read counts data for all samples
all_raw_read_counts_data = {}

# Iterate through each BAM file
for bam_file in bam_files:
    sample_name = os.path.basename(bam_file).split("_")[1].split("-")[0]

    print("Processing:", bam_file)
    print("Sample Name:", sample_name)

    full_bam_file = os.path.join(bam_files_directory, bam_file)

    # Calculate raw read counts for the current BAM file and capture myintrons dictionary
    myintrons = calculate_splicing_efficiency(full_bam_file)

    # Iterate through the myintrons dictionary and format the data
    for start, info in myintrons.items():
        splice_event = f"{start[0]}_{start[1]}"
        unspliced_5p = f"{start[0]}_5p_unspliced"
        unspliced_3p = f"{start[1]}_3p_unspliced"

        # Extract the raw read counts for each event
        raw_spliced_count = info["splice_support"]
        raw_5p_unspliced_count = info["fiveprime_unspliced"] - info["fiveprime_boundary"]
        raw_3p_unspliced_count = info["threeprime_unspliced"] - info["threeprime_boundary"]

        # Update the all_raw_read_counts_data dictionary
        if splice_event not in all_raw_read_counts_data:
            all_raw_read_counts_data[splice_event] = {}
        if unspliced_5p not in all_raw_read_counts_data:
            all_raw_read_counts_data[unspliced_5p] = {}
        if unspliced_3p not in all_raw_read_counts_data:
            all_raw_read_counts_data[unspliced_3p] = {}

        # Store the raw read counts in the dictionary
        all_raw_read_counts_data[splice_event][sample_name] = raw_spliced_count
        all_raw_read_counts_data[unspliced_5p][sample_name] = raw_5p_unspliced_count
        all_raw_read_counts_data[unspliced_3p][sample_name] = raw_3p_unspliced_count

## Create a DataFrame from the raw read counts data
raw_read_counts_df = pd.DataFrame(all_raw_read_counts_data)

# Transpose the DataFrame to switch rows and columns
transposed_raw_df = raw_read_counts_df.T


# Save the transposed DataFrame to a CSV file
output_raw_csv_file = "/homes/biertank/ella/ref_genomes/all_genomes/alignment_output_3/bamfiles/hpv16_bam/all_raw_read_counts_data_transposed.csv"
transposed_raw_df.to_csv(output_raw_csv_file)

completion_message = "Your raw read counts analysis is complete (Yeeeehaw)! You can find the .csv file in your output directory"
emojis = "🎉🎉🎉🎉"
# Print the completion message in green with emojis on both sides
completed_message = f"{emojis} {Fore.GREEN}{completion_message}{Style.RESET_ALL} {emojis}"
print(completed_message)