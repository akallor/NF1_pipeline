#!/usr/bin/env

#conda activate subread_env

featureCounts -T 16 -p -t gene -g gene_id -a ../genome/gencode.v47.annotation.gtf -o gene_counts_exon.txt JH-2-*/*_Aligned.out.bam
