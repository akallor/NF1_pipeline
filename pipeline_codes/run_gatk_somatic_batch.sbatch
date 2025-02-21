#!/bin/bash

#SBATCH --job-name=mutect2_pipeline
#SBATCH --cpus-per-task=8            # Use 8 CPUs per task
#SBATCH --mem=50G                    # Request 32GB RAM
#SBATCH --time=24:00:00               # Max job runtime (adjust as needed)
#SBATCH --output=logs/mutect2_%A_%a.out  # Log file for each job
#SBATCH --error=logs/mutect2_%A_%a.err   # Error log file

source ~/.bashrc

export PATH="~/tools/miniconda3/envs/gatk/bin/java:$PATH"

conda activate gatk

java -version

# Set genome reference directory
gatk="/data/collaborators/aak/NF1_sequencing/gatk-4.5.0.0/gatk"
genome="/data/collaborators/aak/NF1_sequencing/genome"

# Create log directory if it doesn’t exist
mkdir -p logs

# Get list of BAM files
BAM_FILES=($(ls JH-2-*/*.out.sorted.dedup.bam))

# Get the BAM file for this job array index
bamfile=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}

# Extract sample name
sample=$(basename "$bamfile" .out.sorted.dedup.bam)

# Define output file paths
bam_dir=$(dirname "$bamfile")
split_bam="${bam_dir}/${sample}_split.bam"
f1r2_counts="${bam_dir}/${sample}_f1r2.tar.gz"
mutect_raw_vcf="${bam_dir}/${sample}_mutect2_raw.vcf"
mutect_filtered_vcf="${bam_dir}/${sample}_mutect2_filtered.vcf"
orientation_model="${bam_dir}/${sample}_orientation_model.tar.gz"

echo "Processing sample: $sample on node $(hostname)"

# Step 1: Split Reads at Splice Junctions
$gatk --java-options "-Xmx28G" SplitNCigarReads \
    -R "${genome}/hg38.fa" \
    -I "$bamfile" \
    -O "$split_bam"

# Step 2: Collect F1R2 Counts
$gatk --java-options "-Xmx28G" CollectF1R2Counts \
    -R "${genome}/hg38.fa" \
    -I "$split_bam" \
    -O "$f1r2_counts"

# Step 3: Run Mutect2 (Somatic Variant Calling)
$gatk --java-options "-Xmx28G" Mutect2 \
    -R "${genome}/hg38.fa" \
    -I "$split_bam" \
    --germline-resource gnomad.vcf.gz \
    --f1r2-tar-gz "$f1r2_counts" \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --min-base-quality-score 20 \
    -O "$mutect_raw_vcf"

# Step 4: Learn Read Orientation Model
$gatk --java-options "-Xmx28G" LearnReadOrientationModel \
    -I "$f1r2_counts" \
    -O "$orientation_model"

# Step 5: Filter Mutect2 Calls
$gatk --java-options "-Xmx28G" FilterMutectCalls \
    -V "$mutect_raw_vcf" \
    --ob-priors "$orientation_model" \
    -O "$mutect_filtered_vcf"

echo "Finished processing sample: $sample"

