#!/usr/bin/env bash

#Align the reads to the reference indexed genome

./STAR --runThreadN 8 --genomeDir genome/ --readFilesIn 1_S1_L003_R1_001.fastq 1_S1_L003_R2_001.fastq --outFileNamePrefix sample --outSAMtype BAM SortedByCoordinate
