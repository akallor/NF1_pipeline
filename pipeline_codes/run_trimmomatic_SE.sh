#!/usr/bin/env bash
trim_path="/data/collaborators/aak/NF1_sequencing/Trimmomatic-0.39/trimmomatic-0.39.jar"
R1="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data/JH-2-023/GATCTATC-ATGAGGCT_S16_L002_R1_001.fastq"
adapters_path="/data/collaborators/aak/NF1_sequencing/adapters.fa"

# Define output file
out_R1="output_JH-2-023_trimmed.fastq"

java -jar $trim_path SE $R1 \
$out_R1 \
ILLUMINACLIP:$adapters_path:2:30:10 \
SLIDINGWINDOW:4:15 \
MINLEN:36
