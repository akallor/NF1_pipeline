#!/bin/bash



# Iterate over matching BAM files
for bamfile in JH-2-*/*.out.sorted.dedup.bam; do
    samtools index -@ 32 "$bamfile"
done

echo "Indexing complete."

