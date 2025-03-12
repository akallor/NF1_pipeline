#!/bin/bash

# Output file
output_file="reads_mapped_percentage.tsv"
echo -e "file_name\treads_mapped_percentage" > "$output_file"

# Loop through all merged BAM files
for bam in JH-2-*/*_merged.bam; do
    if [[ -f "$bam" ]]; then
        # Extract filename
        filename=$(basename "$bam")

        # Run samtools flagstat and extract the mapped percentage
        mapped_percentage=$(samtools flagstat "$bam" | awk '/mapped \(/ {print $5}' | tr -d '()%')

        # Write results to output file
        echo -e "$filename\t$mapped_percentage" >> "$output_file"
    fi
done

echo "Done! Results saved in $output_file."

