
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
        # Initialize the start_positions dictionary
        start_positions = {}
        formatted_data = []

        for start, info in myintrons.items():
            # Check if splice_support is greater than 1 before including in the analysis
            if info['splice_support'] > 1:
                unspliced_reads = info['fiveprime_unspliced'] + info['threeprime_unspliced']

                # Calculate splicing efficiencies for both 5' and 3' directions
                splicing_efficiency_fiveprime = (info['splice_support'] / (
                        info['splice_support'] + (info['fiveprime_unspliced'] - info['fiveprime_boundary']))) * 100
                splicing_efficiency_threeprime = (info['splice_support'] / (
                        info['splice_support'] + (info['threeprime_unspliced'] - info['threeprime_boundary']))) * 100

                # Update the dictionary entry with rounded splicing efficiencies
                info['splicing_efficiency_fiveprime'] = round(splicing_efficiency_fiveprime, 2)
                info['splicing_efficiency_threeprime'] = round(splicing_efficiency_threeprime, 2)

                # accumulates the total unspliced reads associated with each unique intron start position
                if start[0] not in start_positions:
                    start_positions[start[0]] = unspliced_reads
                else:
                    start_positions[start[0]] += unspliced_reads

                #make df for further datanalysis
                formatted_data.append((
                    start,
                    info['splice_support'],
                    info['fiveprime_unspliced'],
                    info['threeprime_unspliced'],
                    info['fiveprime_boundary'],
                    info['threeprime_boundary'],
                    info['splicing_efficiency_fiveprime'],
                    info['splicing_efficiency_threeprime']
                ))

        # Create a DataFrame from the data
        columns = [
            "Start",
            "Splice Support",
            "5' Unspliced",
            "3' Unspliced",
            "5' Boundary",
            "3' Boundary",
            "5' Splicing Efficiency (%)",
            "3' Splicing Efficiency (%)"
        ]

        pd.options.mode.chained_assignment = None  # Suppress the warning

        # Create a dataframe of the values appended to the formatted_data output
        df = pd.DataFrame(formatted_data, columns=columns)
        # Create a second df copy for the 5prime analysis
        df5 = df

        # Extract the 5' start positions from the "Start" column using .loc
        df5["5' End"] = df["Start"].apply(lambda x: x[0])
        # Drop duplicates based on the "5' End" column, only have unique
        df5 = df5.drop_duplicates(subset="5' End", keep="first")

        # Group by unique start positions for 5' side
        grouped_5prime = df5.groupby("5' End").agg({
            "Splice Support": "sum",
            "5' Unspliced": "sum",
            "5' Boundary": "first"  # Take the first value of 5' Boundary for each group
        }).reset_index()

        # Step 2: Calculate the total unspliced for 5' side
        grouped_5prime["Total Unspliced"] = grouped_5prime["5' Unspliced"] - grouped_5prime["5' Boundary"]

        # Print the compositional analysis for the 5' side
        # print("5' Side Compositional Analysis:")
        # print(grouped_5prime)

        # Filter out rows with the intron value (226, 3357)
        # There is an intron ending at 3360 which has the maximum values, don't need to count both
        df_filtered = df[df["Start"] != (226, 3357)].copy()  # Make a copy of the filtered DataFrame

        # Extract the 3' end positions from the "Start" column using .loc
        df_filtered["3' End"] = df_filtered["Start"].apply(lambda x: x[1])

        # Drop duplicates based on the "3' End" column
        df_filtered = df_filtered.drop_duplicates(subset="3' End", keep="first")

        # Filter out the row with the 3' end position of 3357
        df_filtered = df_filtered[df_filtered["3' End"] != 3357]

        # Step 1: Group by unique 3' end positions for 3' side
        grouped_3prime = df_filtered.groupby("3' End").agg({
            "Splice Support": "sum",
            "3' Boundary": "first",  # Take the first value of 3' Boundary for each group
            "3' Unspliced": "sum"  # Sum the 3' unspliced values for each group
        }).reset_index()

        # Step 2: Calculate the total unspliced for 3' side
        grouped_3prime["Total Unspliced"] = grouped_3prime["3' Unspliced"] - grouped_3prime["3' Boundary"]

        # Print the compositional analysis for the 3' side
        # print("Filtered 3' Side Compositional Analysis:")
        # print(grouped_3prime)

        # Calculate the sum of "Total Unspliced" values for both sides
        # This has to be added to the splice support values that combine to make the denominator for the compositional analysis
        total_unspliced = grouped_5prime["Total Unspliced"].sum() + grouped_3prime["Total Unspliced"].sum()

        # print("Total Unspliced:", total_unspliced)

        # Now calculate the sum of all the splice_support values
        total_splice_support = df['Splice Support'].sum()
        # print("Total Splice Support:", total_splice_support)

        combined_denominator = total_splice_support + total_unspliced
        # print("Combined Denominator:", combined_denominator)

        # Calculate and update the spclearlicing efficiency for each row in the DataFrame
        df['Splicing Efficiency'] = df['Splice Support'] / combined_denominator

        # Select the desired columns from the DataFrame
        desired_columns = [
            'Start', 'Splicing Efficiency', '5\' Unspliced',
            '3\' Unspliced'
        ]
        selected_data = df[desired_columns]
        # # Print the selected data
        # print(selected_data)

        # Create DataFrames for 5' and 3' unspliced percentages
        unique_5prime = grouped_5prime[['5\' End', 'Total Unspliced']]
        unique_3prime = grouped_3prime[['3\' End', 'Total Unspliced']]

        # Calculate the unspliced percentages for 5' and 3' introns
        unique_5prime['5\' Unspliced (%)'] = (unique_5prime['Total Unspliced'] / combined_denominator)
        unique_3prime['3\' Unspliced (%)'] = (unique_3prime['Total Unspliced'] / combined_denominator)

        # Create a column 'Start' with None as the first position for 3' introns
        unique_5prime['Intron'] = unique_5prime['5\' End'].apply(lambda x: (x, "Unspliced 5p"))
        unique_3prime['Intron'] = unique_3prime['3\' End'].apply(lambda x: ("Unspliced 3p", x))

        # Drop the '3\' End' column
        unique_3prime = unique_3prime.drop(columns=['3\' End'])

        # Concatenate the 5' and 3' DataFrames
        combined_unspliced_percentages = pd.concat([unique_5prime, unique_3prime], ignore_index=True)

        # Create a new DataFrame for combined introns and unspliced percentages
        combined_intron_unspliced = pd.DataFrame({
            'Intron': combined_unspliced_percentages['Intron'],
            'Unspliced': combined_unspliced_percentages['5\' Unspliced (%)'].fillna(
                combined_unspliced_percentages['3\' Unspliced (%)'])
        })

        # Print the combined intron and unspliced percentages DataFrame
        # print("Combined Intron and Unspliced Percentages:")
        # print(combined_intron_unspliced)

        # Append the combined_intron_unspliced DataFrame to the bottom of selected_data
        combined_data = pd.concat([selected_data, combined_intron_unspliced], ignore_index=True)

        # Extract the columns needed for the combined data
        combined_data = df[['Start', 'Splicing Efficiency', '5\' Unspliced', '3\' Unspliced']].copy()
        # Append the combined_intron_unspliced DataFrame to the bottom of selected_data
        combined_data = pd.concat([combined_data, combined_intron_unspliced], ignore_index=True)

        # Set NaN values in the 'Start' column with corresponding values from the 'Intron' column, drop 'Intron'.
        combined_data['Start'].fillna(combined_data['Intron'], inplace=True)
        combined_data.drop(columns=['Intron'], inplace=True)
        combined_data['Splicing Efficiency'].fillna(combined_data['Unspliced'], inplace=True)

        combined_data = combined_data[['Start', 'Splicing Efficiency']]
        # Rename the "Splicing Efficiency" column to "Splice Composition"
        combined_data.rename(columns={'Start': "Splice Event"}, inplace=True)
        combined_data.rename(columns={'Splicing Efficiency': sample_name}, inplace=True)

        # # Print the rearranged and combined DataFrame with the renamed column
        # print("Rearranged and Combined DataFrame:")
        # print(combined_data)

        # sum_combined = combined_data[sample_name].sum()
        # print("Total Spliced/Unspliced Composition:", sum_combined)
    return (
        myintrons,
        splice_support,
        non_splice_support,
        total_reads,
        combined_data
    )

# Define the input directory for BAM files
bam_files_directory = "/homes/biertank/ella/ref_genomes/all_genomes/alignment_output_3/bamfiles/hpv16_bam"
# Get a list of BAM files in the specified directory
bam_files = [f for f in os.listdir(bam_files_directory) if f.endswith(".bam")]

# Initialize the dictionary to store splicing data for all samples
all_splicing_data = {}
total_files = len(bam_files)

# Iterate through each BAM file
with tqdm(total=total_files, desc="Processing BAM files", bar_format="{desc} {bar}") as pbar:
    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split("_")[1].split("-")[0]

        print("Processing:", bam_file)
        print("Sample Name:", sample_name)

        full_bam_file = os.path.join(bam_files_directory, bam_file)

        # Calculate splicing efficiency for the current BAM file
        intron_dic, splice_support, non_splice_support, total_reads, combined_data = calculate_splicing_efficiency(full_bam_file)

        # Update the all_splicing_data dictionary
        for _, row in combined_data.iterrows():
            splice_event = row["Splice Event"]
            splicing_efficiency = row[sample_name]  # Adjust the column name

            # Convert the splice_event tuple to a string and replace spaces with underscores
            splice_event_str = "_".join(map(str, splice_event))

            if splice_event_str not in all_splicing_data:
                all_splicing_data[splice_event_str] = {}
            all_splicing_data[splice_event_str][sample_name] = splicing_efficiency

        # Update the loading bar with colorized output
        pbar.set_description(f"{Fore.RED}Processing{Style.RESET_ALL} BAM files")
        pbar.set_postfix(sample=sample_name)
        pbar.update(1)

# Convert the all_splicing_data dictionary to a DataFrame
all_splicing_df = pd.DataFrame(all_splicing_data)
# Transpose the DataFrame to switch rows and columns
transposed_df = all_splicing_df.transpose()
# Rename the index column to "Splice_events"
transposed_df.index.name = "Splice_events"
# Save the transposed DataFrame to a CSV file
output_csv_file = "/homes/biertank/ella/ref_genomes/all_genomes/alignment_output_3/bamfiles/hpv16_bam/all_splicing_data.csv"
transposed_df.to_csv(output_csv_file)

# Define the completion message with science-based emojis
completion_message = "Your initial splice site analysis is complete (Yeeeehaw)! You can find the .csv file in your output directory"
emojis = "🎉🎉🎉🎉"
# Print the completion message in green with emojis on both sides
completed_message = f"{emojis} {Fore.GREEN}{completion_message}{Style.RESET_ALL} {emojis}"
print(completed_message)
