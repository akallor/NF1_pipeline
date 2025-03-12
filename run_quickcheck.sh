#!/bin/bash

# Output file
output_file="quickcheck_results.tsv"

# Write header (optional)
echo -e "File\tStatus" > "$output_file"

# Iterate over matching BAM files
for bamfile in JH-2-*/*.out.sorted.dedup.bam; do
    if samtools quickcheck "$bamfile"; then
        status="ok"
    else
        status="not ok"
    fi
    echo -e "$bamfile\t$status" >> "$output_file"
done

echo "Check complete. Results saved in $output_file"

