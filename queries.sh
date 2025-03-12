#Identity all the clinically pathogenic/likely pathogenic splicing variants from the files

grep '\bNF1\b' *.bed | grep 'splicing' | grep 'pathogenic' > clinvar_NF1_pathogenic_splicing_variants.bed
