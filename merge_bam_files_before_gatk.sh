#!/bin/bash

conda activate gatk
mkdir -p logs
mkdir -p merge_jobs

for sample_dir in JH-2-*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        output_bam="${sample_name}_merged.bam"
        job_script="merge_jobs/merge_${sample_name}.sh"

        cat <<EOF > $job_script
#!/bin/bash
#SBATCH --job-name=merge_${sample_name}
#SBATCH --output=logs/${sample_name}_merge.out
#SBATCH --error=logs/${sample_name}_merge.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G


echo "Merging BAMs for sample: $sample_name"

# Find BAM files safely
bam_files=(\$(find "$sample_dir" -type f -name "*dedup.bam"))

if [ \${#bam_files[@]} -eq 0 ]; then
    echo "No BAMs found for $sample_name — skipping."
    exit 1
fi

# Merge
samtools merge -@ 8 -o ${sample_name}_unsorted.bam "\${bam_files[@]}"

# Sort
samtools sort -@ 8 -o ${sample_name}_merged.bam ${sample_name}_unsorted.bam

# Index
samtools index ${sample_name}_merged.bam

# Cleanup
rm ${sample_name}_unsorted.bam

# Move merged BAM + BAI to sample folder
mv ${sample_name}_merged.bam "$sample_dir/"
mv ${sample_name}_merged.bam.bai "$sample_dir/"

echo "Done merging $sample_name"
EOF

        # Make the job script executable and submit
        chmod +x $job_script
        sbatch $job_script
    fi
done
