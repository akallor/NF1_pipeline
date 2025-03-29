#!/bin/bash
#SBATCH --job-name=annovar_launcher
#SBATCH --output=annovar_launcher.out
#SBATCH --error=annovar_launcher.err
#SBATCH --time=00:10:00  # Short time for job submission script
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

# Set the ANNOVAR directory
ANNOVAR_DIR="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data/annovar"
DB_DIR="humandb/"

# Loop through all JH-2-* directories
for sample_dir in ../JH-2-*; do
    sample_name=$(basename "$sample_dir")

    # Check if it's a directory
    if [[ -d "$sample_dir" ]]; then
        echo "Submitting job for $sample_name"
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=annovar_$sample_name
#SBATCH --output=annovar_$sample_name.out
#SBATCH --error=annovar_$sample_name.err
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=10

module load perl  # Load Perl module if needed

for vcf_file in "$sample_dir"/*filtered.vcf; do
    if [[ -f "\$vcf_file" && "\$vcf_file" != *germline* ]]; then
        echo "Processing VCF file: \$vcf_file"

        # Convert VCF to ANNOVAR input format
        perl "$ANNOVAR_DIR/convert2annovar.pl" -format vcf4 "\$vcf_file" > "\$vcf_file.avinput"

        # Annotate using ANNOVAR
        perl "$ANNOVAR_DIR/table_annovar.pl" "\$vcf_file.avinput" "$DB_DIR" \\
            -buildver hg38 -out "\$vcf_file.annotated" -remove \\
            -protocol refGene,avsnp150,clinvar \\
            -operation g,f,f -nastring .

        echo "Annotation complete for: \$vcf_file"
    fi
done
EOF
    fi
done
