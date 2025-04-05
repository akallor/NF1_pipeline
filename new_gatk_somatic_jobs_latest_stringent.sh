#!/bin/bash
# run_gatk_pipeline.sh
# Load environment and set Java 17
source ~/.bashrc
export JAVA_HOME=~/tools/miniconda3/envs/gatk
export PATH=$JAVA_HOME/bin:$PATH
# Set directories and executables
gatk="/data/collaborators/aak/NF1_sequencing/gatk-4.5.0.0/gatk"
bcftools="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data/bcftools-1.21/bcftools/bcftools"
genome="/data/collaborators/aak/NF1_sequencing/genome"
reference="${genome}/hg38.fa"  # adjust as needed
known_site="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data/known_site"
# Known sites for BQSR (typically dbSNP, 1000 Genomes, Mills indels)
#known_sites_1="${known_site}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
#known_sites_2="${known_site}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
known_sites_dbsnp="${known_site}/updated_known_sites.vcf"

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
#SBATCH --output=${sample_dir}/${bam_base}_GATK.out
#SBATCH --error=${sample_dir}/${bam_base}_GATK.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
# Load environment in the job
source ~/.bashrc
export JAVA_HOME=~/tools/miniconda3/envs/gatk
export PATH=\$JAVA_HOME/bin:\$PATH
echo "Processing BAM file: ${bam_file} in sample directory: ${sample_dir}"

# 0. Base Quality Score Recalibration (BQSR)
# First pass: BaseRecalibrator to create recalibration table
echo "Running BaseRecalibrator (First Pass)..."
$gatk BaseRecalibrator \
    -R ${reference} \
    -I ${bam_file} \
    --known-sites ${known_sites_dbsnp} \
    -O ${sample_dir}/${bam_base}.recal_data.table

# Apply Base Quality Score Recalibration
echo "Applying Base Quality Score Recalibration..."
$gatk ApplyBQSR \
    -R ${reference} \
    -I ${bam_file} \
    --bqsr-recal-file ${sample_dir}/${bam_base}.recal_data.table \
    -O ${sample_dir}/${bam_base}.recalibrated.bam

# Optional: Second pass BaseRecalibrator to assess improvement
echo "Running BaseRecalibrator (Second Pass)..."
$gatk BaseRecalibrator \
    -R ${reference} \
    -I ${sample_dir}/${bam_base}.recalibrated.bam \
    --known-sites ${known_sites_dbsnp} \
    -O ${sample_dir}/${bam_base}.post_recal_data.table

# Generate BQSR report (optional, but recommended)
echo "Generating BQSR Report..."
$gatk AnalyzeCovariates \
    -before ${sample_dir}/${bam_base}.recal_data.table \
    -after ${sample_dir}/${bam_base}.post_recal_data.table \
    -plots ${sample_dir}/${bam_base}.recalibration_plot.pdf

# 1. SplitNCigarReads (use recalibrated BAM)
echo "Running SplitNCigarReads..."
$gatk SplitNCigarReads \
    --reference ${reference} \
    --input ${sample_dir}/${bam_base}.recalibrated.bam \
    --output ${sample_dir}/${bam_base}.split.bam

# 2. Collect F1R2 Counts
echo "Collecting F1R2 Counts..."
$gatk CollectF1R2Counts \
    -I ${sample_dir}/${bam_base}.split.bam -R ${reference}\
    -O ${sample_dir}/${bam_base}.f1r2.tar.gz

# 3. Mutect2: Somatic variant calling
echo "Running Mutect2..."
$gatk Mutect2 \
    -R ${reference} \
    -I ${sample_dir}/${bam_base}.split.bam \
    --f1r2-tar-gz ${sample_dir}/${bam_base}.f1r2.tar.gz \
    -O ${sample_dir}/${bam_base}.unfiltered.vcf

# 4. Learn Read Orientation Model
echo "Learning read orientation model..."
$gatk LearnReadOrientationModel \
    -I ${sample_dir}/${bam_base}.f1r2.tar.gz \
    -O ${sample_dir}/${bam_base}.read-orientation-model.tar.gz

# 5. Filter Mutect2 Calls (Stricter Filters)
echo "Filtering Mutect2 calls with stricter thresholds..."
$gatk FilterMutectCalls \
    -V ${sample_dir}/${bam_base}.unfiltered.vcf \
    --ob-priors ${sample_dir}/${bam_base}.read-orientation-model.tar.gz \
    --min-median-base-quality 20 \
    --max-alt-allele-count 3 \
    -O ${sample_dir}/${bam_base}.filtered.vcf -R ${reference}

# 6. Additional VariantFiltration (Hard Filters)
echo "Applying hard filters to remove false positives..."
$gatk VariantFiltration \
    -R ${reference} \
    -V ${sample_dir}/${bam_base}.filtered.vcf \
    --filter-expression "DP < 20" --filter-name "LowDepth" \
    -O ${sample_dir}/${bam_base}.hardfiltered.hc.vcf

# 7. Keep only "PASS" variants
echo "Extracting only PASS variants..."
$bcftools view -f "PASS" ${sample_dir}/${bam_base}.hardfiltered.vcf > ${sample_dir}/${bam_base}.final.hc.vcf

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
