#!/usr/bin bash
#Correct the BQSR recalibration error due to mismatched headers

#Run the following commands to inspect contig names:

samtools view -H your_input.bam | grep "@SQ"

#Reference genome
grep ">" reference.fasta | cut -d " " -f1

#If the BAM file lacks "chr" prefixes but the reference has them, use:

samtools reheader -c 'sed s/^@SQ\tSN:\([0-9XY]\)/@SQ\tSN:chr\1/' input.bam > output.bam

#If the reference lacks "chr" prefixes but the BAM file has them, reformat the reference using:

sed -i 's/^>chr/>/' reference.fasta

#After fixing, index the new BAM:
samtools index output.bam

#Modify the chromosome names in your VCF file

#Create a mapping file (rename.txt)

printf "1\tchr1\n2\tchr2\n3\tchr3\n4\tchr4\n5\tchr5\n6\tchr6\n7\tchr7\n8\tchr8\n9\tchr9\n10\tchr10\n11\tchr11\n12\tchr12\n13\tchr13\n14\tchr14\n15\tchr15\n16\tchr16\n17\tchr17\n18\tchr18\n19\tchr19\n20\tchr20\n21\tchr21\n22\tchr22\nX\tchrX\nY\tchrY\nMT\tchrM\n" > rename.txt

#Edit the VCF header to explicitly list the contigs (chromosomes) before running bcftools annotate

bcftools view -h common_all_20180418.vcf > header.vcf

#Append the following lines to header.vcf

##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=X>
##contig=<ID=Y>
##contig=<ID=MT>

#Update your VCF with the corrected header

bcftools reheader -h header.vcf -o fixed_known_sites.vcf common_all_20180418.vcf

#Rename the chromosomes using your rename.txt file

bcftools annotate --rename-chrs rename.txt -o updated_known_sites.vcf fixed_known_sites.vcf

#Remember: Before running GATK, always have indexed and dictionary files



