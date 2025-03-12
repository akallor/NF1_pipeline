#!/bin/bash

# Define paths
SCHEMA="schema2.tsv"
BBMAP="/data/collaborators/aak/NF1_sequencing/bbmap"  # Update with correct BBMap path
OUTPUT_BASE="/data/collaborators/aak/NF1_sequencing/John_Hopkins_data"
LOG_DIR="${OUTPUT_BASE}/logs"

source ~/.bashrc

mkdir -p "$LOG_DIR"  # Ensure log directory exists

# Read schema and group files by patient and lane/index
declare -A file_groups

while IFS=$'\t' read -r patient file; do
    mkdir -p "${OUTPUT_BASE}/${patient}"  # Ensure patient directory exists
    
    # Extract unique key for grouping (e.g., lane, index, sample ID)
    key=$(echo "$file" | sed -E 's/(_R[12]_trimmed.fastq)//')
    
    # Store files in an array under the patient+key identifier
    file_groups["${patient},${key}"]+="$file "
done < <(tail -n +2 "$SCHEMA")  # Skip the header line

# Submit a SLURM job for each unique pair
for group in "${!file_groups[@]}"; do
    patient=$(echo "$group" | cut -d',' -f1)
    key=$(echo "$group" | cut -d',' -f2)

    # Extract R1 and R2
    fq1=$(echo "${file_groups[$group]}" | tr ' ' '\n' | grep '_R1_' | head -n1)
    fq2=$(echo "${file_groups[$group]}" | tr ' ' '\n' | grep '_R2_' | head -n1)

    if [[ -f "${OUTPUT_BASE}/${patient}/$fq1" && -f "${OUTPUT_BASE}/${patient}/$fq2" ]]; then
        out1="${OUTPUT_BASE}/${patient}/repair_${fq1}"
        out2="${OUTPUT_BASE}/${patient}/repair_${fq2}"

        # Create SLURM job script
        job_script="${OUTPUT_BASE}/${patient}/run_bbmap_repair_${patient}_${key}.sh"
        cat <<EOF > "$job_script"
#!/bin/bash
#SBATCH --job-name=BBMAP_${patient}_${key}
#SBATCH --output=${LOG_DIR}/bbmap_${patient}_${key}_%j.out
#SBATCH --error=${LOG_DIR}/bbmap_${patient}_${key}_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

echo "Running BBMap repair for $patient ($key)..."

$BBMAP/repair.sh in1="${OUTPUT_BASE}/${patient}/$fq1" in2="${OUTPUT_BASE}/${patient}/$fq2" out1="$out1" out2="$out2" repair

if [[ \$? -eq 0 ]]; then
    echo "Repair completed for $patient ($key). Output saved in ${OUTPUT_BASE}/${patient}/"
else
    echo "Error running BBMap repair for $patient ($key)!"
fi
EOF

        # Make the script executable
        chmod +x "$job_script"

        # Submit the job to SLURM
        sbatch "$job_script"
    else
        echo "Warning: Missing paired files for $patient ($key), skipping submission."
    fi
done

