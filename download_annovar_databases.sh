
#!/usr/bin bash

#Download the gnomad genome database

./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/


#Download the gnomad exome database

./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_exome humandb/
