#!/bin/bash


BASE_DIR="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data"
samtools="/data/collaborators/aak/NF1_sequencing/samtools"
picard="/data/collaborators/aak/NF1_sequencing/picard.jar"

# Iterate through each directory matching 'JH-2-*'
for dir in "$BASE_DIR"/JH-2-*; do
    [ -d "$dir" ] || continue  # Skip if not a directory

    echo "Processing directory: $dir"

    # Find all BAM files matching '*_Aligned.out.bam'
    for bam in "$dir"/*_Aligned.out.bam; do
        [ -f "$bam" ] || continue  # Skip if no matching BAM files are found

        # Define job name based on BAM file
        job_name=$(basename "$bam" .bam)

        # Submit a separate SLURM job for each BAM file
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$dir/${job_name}.log
#SBATCH --error=$dir/${job_name}.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

echo "Processing BAM file: $bam"

# Sort BAM file
sorted_bam="${bam%.bam}.sorted.out.bam"
$samtools sort -o "\$sorted_bam" "$bam"

# Index the sorted BAM file
$samtools index "\$sorted_bam"

# Run Picard to remove duplicates
dedup_bam="${sorted_bam%.sorted.out.bam}.dedup.bam"
metrics_file="${sorted_bam%.sorted.out.bam}.duplication_metrics.txt"

java -jar $picard MarkDuplicates \
    -I "\$sorted_bam" \
    -O "\$dedup_bam" \
    -M "\$metrics_file" \
    -REMOVE_DUPLICATES true \
    -CREATE_INDEX true

echo "Finished processing: $bam"
EOF

    done
done

echo "All SLURM jobs submitted!"

