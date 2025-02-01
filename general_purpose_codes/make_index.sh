#!/usr/bin/env bash

#Create the index of the reference genome

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles genome/hg38.fa
