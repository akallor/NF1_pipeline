#!/usr/bin/env


#Identifying Circular RNAs with FindCirc

find_circ.py -r reference.fasta -1 read1.fastq -2 read2.fastq -o output_prefix

#Annotating Non-Coding RNAs with HOMER 

#Install HOMER: Download and install HOMER from its official website 

#Prepare Your Annotation File: Ensure you have a GTF file with annotations for non-coding RNAs 

#Run HOMER Annotation: Use the following command to annotate your reads 

annotatePeaks.pl aligned.bam gtf_annotation.gtf gtf_annotation > annotated_peaks.txt 
