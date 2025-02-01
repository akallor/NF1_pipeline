#!/usr/bin/env bash

#Perform QC of all the patient data files in the directory

fastqc JH-*

#Perform adapter trimming using trimmomoatic

#Perform indexing of the genome

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles genome/hg38.fa 

#Perform read alignment

./STAR --runThreadN 8 --genomeDir genome/ --readFilesIn 1_S1_L003_R1_001.fastq 1_S1_L003_R2_001.fastq --outFileNamePrefix sample --outSAMtype BAM SortedByCoordinate

#Perform transcript assembly

./stringtie sampleAligned.sortedByCoord.out.bam -o sample_transcripts.gtf -p 8 

#Normalize and quantify the reads
