#!/usr/bin python

import pandas as pd
import glob
import os

# Find all files matching the pattern
file_pattern = "JH-2-*/*.final.vcf.annotated.hg38_multianno.txt"
all_files = glob.glob(file_pattern)

# Check if files were found
if not all_files:
    print(f"No files found matching pattern: {file_pattern}")
    exit(1)

print(f"Found {len(all_files)} files matching the pattern.")

# Get unique file names (in case of duplicates)
unique_files = []
seen_files = set()

for file_path in all_files:
    # Get just the file name without the path
    file_name = os.path.basename(file_path)
    if file_name not in seen_files:
        seen_files.add(file_name)
        unique_files.append(file_path)

print(f"After removing duplicates: {len(unique_files)} unique files to process.")

# Initialize an empty list to store dataframes
dfs = []

# Process each file
for i, file_path in enumerate(unique_files):
    try:
        # Read the file (assuming it's tab-separated)
        df = pd.read_csv(file_path, sep='\t')

        # Add file information for debugging if needed
        # df['source_file'] = os.path.basename(file_path)

        # Append to our list
        dfs.append(df)

        # Print progress
        if (i + 1) % 10 == 0 or i + 1 == len(unique_files):
            print(f"Processed {i + 1}/{len(unique_files)} files...")

    except Exception as e:
        print(f"Error processing {file_path}: {e}")

# Concatenate all dataframes
if dfs:
    result = pd.concat(dfs, ignore_index=True)

    # Save to the output file
    output_file = "complete_annotation_latest.tsv"
    result.to_csv(output_file, sep='\t', index=False)

    print(f"Successfully concatenated {len(dfs)} files into {output_file}")
    print(f"Total rows in output file: {len(result)}")
else:
    print("No valid data was found to concatenate.")
