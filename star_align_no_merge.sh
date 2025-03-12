#!/bin/bash

# Define paths
GENOME_DIR="/data/collaborators/aak/NF1_sequencing/genome"
BASE_DIR="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data/testing"
STAR="/data/collaborators/aak/NF1_sequencing/STAR"
samtools="/data/collaborators/aak/NF1_sequencing/samtools-1.21/samtools"

#SBATCH --job-name=STAR_Alignment
#SBATCH --output=alignment_output.log
#SBATCH --error=alignment_error.log
#SBATCH --time=6:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1

for subdir in JH-2-*; do
    if [ -d "$subdir" ]; then
        cd "$subdir"
        echo "Processing directory: $subdir"  # Debugging line

        # Collect job IDs for STAR alignments
        STAR_JOBS=()
        for file in *_1.fastq; do
            echo "Found file: $file"  # Debugging line

            PREFIX=$(basename "$file" | sed 's/_1.fastq//')
            R2_FILE="${file/_1/_2}"
            echo "Looking for pair: $R2_FILE"  # Debugging line

            # Ensure R2 file exists before submitting job
            if [ ! -f "$R2_FILE" ]; then
                echo "Error: Missing pair file for $file"
                continue
            fi

            STAR_JOB_ID=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=STAR_${PREFIX}
#SBATCH --output=${PREFIX}_STAR.out
#SBATCH --error=${PREFIX}_STAR.err
#SBATCH --time=6:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1

$STAR --runThreadN 8 --genomeDir "$GENOME_DIR" --readFilesIn "$file" "$R2_FILE" --outFileNamePrefix "${PREFIX}_" --outSAMtype BAM Unsorted
EOF
            )
            STAR_JOBS+=("$STAR_JOB_ID")
        done

        # Ensure all STAR jobs complete before merging
        if [ ${#STAR_JOBS[@]} -eq 0 ]; then
            echo "No STAR jobs were submitted for $subdir, skipping merge."
            cd ..  # Go back to the base directory
            continue
        fi

        JOB_IDS=$(IFS=,; echo "${STAR_JOBS[*]}")
        echo "Waiting for STAR jobs to finish before merging in $subdir"

        # Submit samtools merge job with explicit dependency
        sbatch --dependency=afterok:$JOB_IDS <<EOF
#!/bin/bash
#SBATCH --job-name=MergeBAMs_${subdir}
#SBATCH --output=merge_bams.out
#SBATCH --error=merge_bams.err
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Wait until all BAM files exist before merging
for ((i=0; i<30; i++)); do
    if ls *_Aligned.out.bam &>/dev/null; then
	echo "BAM files detected for merging."
        break
    fi
    sleep 10
    echo "Waiting for BAM files in $subdir..."
done

if ! ls *_Aligned.out.bam &>/dev/null; then
    echo "Error: No BAM files found in $subdir after waiting. Exiting merge."
    exit 1
fi


# Sort and index final BAM
$samtools sort -@ 4 -o "${subdir}_sorted.bam" "${subdir}.bam"
$samtools index "${subdir}_sorted.bam"
EOF

        cd ..  # Go back to the base directory
    fi
done

