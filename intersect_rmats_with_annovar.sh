#!/usr/bin bash





#Get all the annotated vcf files from somatic filtered annotated files by annovar

awk 'FNR>1||NR==1' ../JH-2-*/*.filtered.vcf.annotated.hg38_multianno.txt > all_jh2.txt

#Relaxed FDR filter for rmats output (to get NF1 splicing events) (5% instead of 1%) 

awk -F'\t' 'NR==1 || ($20 <= 5) {print $4, $6, $7, $19, $20}' OFS='\t' */*.MATS.JC.txt > all_conditions_5perc.bed

#Re-intersect with reduced filter

bedtools intersect -a all_jh2.txt -b all_conditions_5perc.bed -wa -wb > clinvar_overlaps_5perc.bed
