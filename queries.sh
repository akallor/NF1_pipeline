#Identity all the clinically pathogenic/likely pathogenic splicing variants from the files

#Updated query

grep '\bNF1\b' *.bed | grep 'splicing' | egrep 'Pathogenic|pathogenic' | awk '{print $13, $14, $15,$6,$7,$12,$8, $9, $2, $3}' OFS="\t" > clinvar_NF1_pathogenic.bed

#Extract nucleotide sequences from the bed file to prepare for Transdecoder

bedtools getfasta -fi ../../genome/hg38.fa -bed clinvar_NF1_pathogenic.bed -fo extracted_NF1_clinvar_sequences.fa

#Get all splice-site variants from all the filtered and annotated txt files in all samples

grep '\bNF1\b' JH-2-*/*.hg38_multianno.txt | grep '\bsplicing\b' | egrep 'Pathogenic|pathogenic' | awk -F':' '{print $2,$6}
