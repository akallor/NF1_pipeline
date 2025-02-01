#!/usr/bin/env bash

#Run featurecounts to check how many reads map to genomic features

featureCounts -a annotation.gtf -o counts.txt -t exon -g gene -p aligned.bam
