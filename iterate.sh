#!/bin/bash

# Path to the Python script
python_script="/homes/biertank/ella/ref_genomes/all_genomes/splice/final_splicesites.py"

# Output directory
output_dir="/homes/biertank/ella/ref_genomes/all_genomes/splice/csv_splice"

# Loop through BAM files
for bam_file in /homes/biertank/ella/ref_genomes/all_genomes/hpv_alignment/HPV_16/clustering/*.bam; do
    # Run the Python script
    python "$python_script" "$bam_file" "$output_dir" 
done

