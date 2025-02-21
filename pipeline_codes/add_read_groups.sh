#!/bin/bash

BASE_DIR="/data/collaborators/aak/NF1_sequencing"
picard="${BASE_DIR}/picard.jar"


# SLURM parameters
MEMORY="16G"
CPUS="4"
TIME="1:00:00"

for dir in JH-2-*; do
    if [[ -d "$dir" ]]; then
        echo "Checking directory: $dir"

        for bamfile in "$dir"/*.out.sorted.out.bam; do
            if [[ -f "$bamfile" ]]; then
                rg_bam="${bamfile/.out.sorted.out.bam/.out.sorted.rg.bam}"
                job_name="add_rg_${dir}_$(basename "$bamfile" .out.sorted.out.bam)"

                # Submit SLURM job to add read groups
                sbatch --mem=$MEMORY --cpus-per-task=$CPUS --time=$TIME \
                    --job-name="$job_name" --output="${bamfile%.bam}_add_rg.log" --wrap="
                    java -jar $picard AddOrReplaceReadGroups \
                        -I '$bamfile' \
                        -O '$rg_bam' \
                        -RGID '$(basename "$bamfile")' \
                        -RGLB 'Library1' \
                        -RGPL 'ILLUMINA' \
                        -RGPU 'Unit1' \
                        -RGSM 'Sample1' \
                        -SORT_ORDER coordinate
                "
                echo "Submitted AddOrReplaceReadGroups job for: $bamfile"
            fi
        done
    fi
done

