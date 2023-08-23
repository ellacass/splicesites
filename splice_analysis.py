import pysam
from sys import argv
from tabulate import tabulate
import csv
import os
import pandas as pd


def calculate_splicing_efficiency(bam_file):
    # Initialize variables to keep track of reads and coverage
    total_reads = 0

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

            ###THIS CALCULATES THE SPLICE SUPPORT READS
            ###CONFIRMED BY BAM FILE EXAMINATION THAT THIS IS CALCULATING CORRECTLY
            ###USES THE SPLUCED LENGTH VAR WHICH CONTAINS READS WITH CIGAR STRINGS "N"
            if spliced_length > 0:
                running_over = [t for t in intron_regions if
                                read.reference_start < t[0] + 3 and read.reference_end > t[1] - 3]
                if len(running_over) > 0:
                    for t in running_over:
                        intron_length = t[1] - t[0]
                        if spliced_length <= intron_length <= spliced_length:
                            myintrons[t]["splice_support"] += 1
                    splice_support += 1

            ###THIS CALCULATES THE FIVE AND THREE PRIME UNSPLICED READS AND BOUNDARY READS
            ###CANNOT TELL IF BOUNDARY READS ARE SPLICED/UNSPLICED, THEREFORE MUST SUBTRACTED FROM UNSPLICED (?)
            if spliced_length == 0:
                fiveprime_unspliced = [t for t in intron_regions if
                                       read.reference_start < t[0] and read.reference_end > t[0]
                                       and read.reference_end < t[1] and
                                       "N" not in read.cigarstring]

                threeprime_unspliced = [t for t in intron_regions if
                                        read.reference_start < t[1] and read.reference_end > t[1] and
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
                                  read.reference_start >= t[1] - 5 and
                                  read.reference_start <= t[1] -1 and
                                  read.reference_end >= t[1] and
                                  "N" not in read.cigarstring]
                if len(three_boundary) > 0:
                    for t in three_boundary:
                        myintrons[t]["threeprime_boundary"] += 1

    return (
        myintrons,
        splice_support,
        non_splice_support,
        total_reads,
    )

bam_file = argv[1]  # Example BAM file name
output_dir = "/homes/biertank/ella/ref_genomes/all_genomes/splice/csv_splice"
intron_dic, splice_support, non_splice_support, total_reads = calculate_splicing_efficiency(bam_file)
sample_name = os.path.basename(bam_file).split("_")[1].split("-")[0]





# Print the input BAM file name
# output_csv_file = f"{sample_name}.csv"  # Define the output CSV file name based on the sample name
print(f"Processing BAM file: {bam_file}")

# Initialize the start_positions dictionary
start_positions = {}
formatted_data=[]
data=[]

# # Loop through intron_dic and update the values
for start, info in intron_dic.items():
    unspliced_reads = info['fiveprime_unspliced'] + info['threeprime_unspliced']
    uncat_reads = info['fiveprime_boundary'] + info['threeprime_boundary']

    # Calculate splicing efficiencies for both 5' and 3' directions
    splicing_efficiency_fiveprime = (info['splice_support'] / (
            info['splice_support'] + (info['fiveprime_unspliced'] - info['fiveprime_boundary']))) * 100
    splicing_efficiency_threeprime = (info['splice_support'] / (
            info['splice_support'] + (info['threeprime_unspliced'] - info['threeprime_boundary']))) * 100

    # Update the dictionary entry with rounded splicing efficiencies
    info['splicing_efficiency_fiveprime'] = round(splicing_efficiency_fiveprime, 2)
    info['splicing_efficiency_threeprime'] = round(splicing_efficiency_threeprime, 2)

    if start[0] not in start_positions:
        start_positions[start[0]] = unspliced_reads
    else:
        start_positions[start[0]] += unspliced_reads
    #
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
#
#
splicing_efficiencies={}
#
# # Open the CSV file for writing
output_csv_file = os.path.join(output_dir, f"{sample_name}.csv")
print("Constructed output_csv_file:", output_csv_file)
print("Sample name:", sample_name)
csv_file = open(output_csv_file, mode='w', newline='')
csv_writer = csv.writer(csv_file)
#
# # Write the header row
csv_writer.writerow(["Sample", "Intron", "5' Splice Eff", "3' Splice Eff"])
#
# # Write the data rows
for intron, efficiencies in splicing_efficiencies.items():
    csv_writer.writerow([sample_name, f"{intron[0]}-{intron[1]}", f"{efficiencies['5prime']}%", f"{efficiencies['3prime']}%"])

# Close the CSV file
csv_file.close()

print(f"CSV file '{output_csv_file}' created.")
# # Print formatted data using tabulate
headers = [
    "Start",
    "Splice Support",
    "5' Unspliced",
    "3' Unspliced",
    "5' Boundary",
    "3' Boundary",
    "5' Splicing Efficiency (%)",
    "3' Splicing Efficiency (%)"
]
#
# print(tabulate(formatted_data, headers=headers, tablefmt="pretty"))
# print(formatted_data)


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


df = pd.DataFrame(formatted_data, columns=columns)
df5=df

# Extract the 3' end positions from the "Start" column using .loc
df5["5' End"] = df["Start"].apply(lambda x: x[0])
# Drop duplicates based on the "3' End" column
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

print("Total Unspliced:", total_unspliced)


# Now calculate the sum of all the splice_support values
total_splice_support = df['Splice Support'].sum()
print("Total Splice Support:", total_splice_support)

combined_denominator = total_splice_support + total_unspliced
print("Combined Denominator:", combined_denominator)

# Calculate and update the splicing efficiency for each row in the DataFrame
df['Splicing Efficiency'] = df['Splice Support'] / combined_denominator

# # Print the updated DataFrame
# print("Updated DataFrame:")
# print(df)

# Select the desired columns from the DataFrame
desired_columns = [
    'Start', 'Splice Support', '5\' Unspliced', '3\' Unspliced',
    '5\' Boundary', '3\' Boundary', 'Splicing Efficiency'
]
selected_data = df[desired_columns]

# Multiply the splicing efficiency values by 100 to represent them as percentages
selected_data['Splicing Efficiency'] = selected_data['Splicing Efficiency'] * 100

# Print the selected data
print(selected_data)

# Sum all values in the "Splicing Efficiency" column
total_splicing_efficiency = selected_data['Splicing Efficiency'].sum()
# Print the total splicing efficiency
print("Total Splicing Efficiency:", total_splicing_efficiency)