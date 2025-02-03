#!?usr/bin/env bash

trim_path="/data/collaborators/aak/NF1_sequencing/Trimmomatic-0.39/trimmomatic-0.39.jar"

R1="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data/JH-2-023/GATCTATC-ATGAGGCT_S16_L002_R1_001.fastq"

#$R2="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data/JH-2-001/CCDYYANXX_8_ATTACTCG-TATAGCCT_2.fastq"


adapters_path="/data/collaborators/aak/NF1_sequencing/adapters.fa"

#java -jar $trim_path PE $R1 $R2 output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz ILLUMINACLIP:overrep.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

java -jar $trim_path PE $R1 $R2 \
$out_R1_paired $out_R1_unpaired \
$out_R2_paired $out_R2_unpaired \
ILLUMINACLIP:$adapters_path:2:30:10 \ 
SLIDINGWINDOW:4:15 \
MINLEN:36 
