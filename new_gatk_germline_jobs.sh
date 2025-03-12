#!/bin/bash
# run_gatk_pipeline.sh

# Load environment and set Java 17
source ~/.bashrc
export JAVA_HOME=~/tools/miniconda3/envs/gatk
export PATH=$JAVA_HOME/bin:$PATH

# Set directories and executables
gatk="/data/collaborators/aak/NF1_sequencing/gatk-4.5.0.0/gatk"
genome="/data/collaborators/aak/NF1_sequencing/genome"
reference="${genome}/hg38.fa"  # adjust as needed

# Loop over directories starting with "JH-2-"
for sample_dir in JH-2-*; do
    if [ -d "$sample_dir" ]; then
        # Find all dedup bam files in the directory
        for bam_file in "$sample_dir"/*.dedup.bam; do
            if [ -f "$bam_file" ]; then
                # Extract base name for output files
                bam_base=$(basename "$bam_file" .dedup.bam)
                echo "Submitting job for BAM file: $bam_file in sample directory: $sample_dir"

                # Submit a SLURM job using a heredoc
                sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${bam_base}_GATK
#SBATCH --output=${sample_dir}/${bam_base}_GATK_haplo.out
#SBATCH --error=${sample_dir}/${bam_base}_GATK_haplo.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load environment in the job
source ~/.bashrc
export JAVA_HOME=~/tools/miniconda3/envs/gatk
export PATH=\$JAVA_HOME/bin:\$PATH

echo "Processing BAM file: ${bam_file} in sample directory: ${sample_dir}"

# 1. SplitNCigarReads
echo "Running SplitNCigarReads..."
$gatk SplitNCigarReads \
    -R ${reference} \
    -I ${bam_file} \
    -O ${sample_dir}/${bam_base}.split.bam

# 2. Run HaplotypeCaller (Germline Variant Calling)
echo "Running HaplotypeCaller..."
$gatk HaplotypeCaller \
    -R ${reference} \
    -I ${sample_dir}/${bam_base}.split.bam \
    -O ${sample_dir}/${bam_base}.germline_raw.vcf

# 3. Filter Germline Calls
echo "Filtering Germline calls..."
$gatk VariantFiltration \
    -V ${sample_dir}/${bam_base}.germline_raw.vcf \
    -O ${sample_dir}/${bam_base}.germline_filtered.vcf \
    --filter-name "QD_filter" --filter-expression "QD < 2.0" \
    --filter-name "FS_filter" --filter-expression "FS > 60.0" \
    --filter-name "MQ_filter" --filter-expression "MQ < 40.0"

echo "Job for BAM file ${bam_file} in sample directory ${sample_dir} completed."
EOF
            else
                echo "No dedup BAM files found in $sample_dir."
            fi
        done
    else
        echo "Directory $sample_dir does not exist."
    fi
done

