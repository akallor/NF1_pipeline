#!?usr/bin/env bash

#Assemble reads into transcripts using Stringtie

./stringtie sampleAligned.sortedByCoord.out.bam -o sample_transcripts.gtf -p 8
