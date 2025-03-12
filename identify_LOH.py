#!/usr/bin/env

#Identifying loss of heterozygosity from vcf files

#Pre-processing step:
#Identify ref, alt and %af from the vcf file (using one of the JH-2-038 vcf files as example)

#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' tumor.vcf > jh2038_af.bed


import pandas as pd

# Read the file correctly (assuming no header in the input file)
jh = pd.read_csv("jh2038_af.bed", sep="\t", header=None, names=['chr', 'pos', 'ref', 'alt', 'perc_alt'])

# Split allele depth into two columns: Reference and Alternate read counts
jh[['perc_alt_ref', 'perc_alt_alt']] = jh['perc_alt'].str.split(',', expand=True)

# Convert to numeric (handle potential errors like missing values)
jh['perc_alt_ref'] = pd.to_numeric(jh['perc_alt_ref'], errors='coerce').fillna(0).astype(int)
jh['perc_alt_alt'] = pd.to_numeric(jh['perc_alt_alt'], errors='coerce').fillna(0).astype(int)

# Compute total reads
jh['total'] = jh['perc_alt_ref'] + jh['perc_alt_alt']

# Compute Variant Allele Frequency (VAF), avoiding division by zero
jh['vaf'] = jh['perc_alt_alt'] / jh['total']
jh['vaf'] = jh['vaf'].fillna(0)  # Replace NaN with 0 if total reads are 0

# Display first few rows
print(jh.head())

