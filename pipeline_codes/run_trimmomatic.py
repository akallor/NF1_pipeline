#!/usr/bin/env python

import pandas as pd
import subprocess
import os
import argparse
from pathlib import Path
from glob import glob

def run_fastqc(input_file, output_dir):
    """Run FastQC on a single file"""
    try:
        cmd = [
            'fastqc',
            '--outdir', output_dir,
            input_file
        ]
        subprocess.run(cmd, check=True)
        print(f"Successfully ran FastQC on {input_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running FastQC on {input_file}: {e}")

def run_trimmomatic(input_fastq, adapter_type, adapters_fa, trimmomatic_path):
    # Construct input/output filenames
    base_name = Path(input_fastq).stem
    if base_name.endswith('_R1'):
        base_name = base_name[:-3]
    elif base_name.endswith('_1'):
        base_name = base_name[:-2]
    
    input_r1 = input_fastq
    input_r2 = input_fastq.replace('_R1', '_R2').replace('_1.', '_2.')
    
    output_dir = 'trimmed'
    os.makedirs(output_dir, exist_ok=True)
    
    output_r1_paired = f"{output_dir}/{base_name}_R1_paired.fastq.gz"
    output_r1_unpaired = f"{output_dir}/{base_name}_R1_unpaired.fastq.gz"
    output_r2_paired = f"{output_dir}/{base_name}_R2_paired.fastq.gz"
    output_r2_unpaired = f"{output_dir}/{base_name}_R2_unpaired.fastq.gz"
    
    # Determine which adapter sequence to use
    adapter_arg = ""
    if adapter_type == "Illumina Universal Adapter":
        adapter_arg = f"ILLUMINACLIP:{adapters_fa}:2:30:10:2:keepboth"
    elif adapter_type == "TruSeq Adapter":
        adapter_arg = f"ILLUMINACLIP:{adapters_fa}:2:30:10:2:keepboth"
    
    # Construct Trimmomatic command with java -jar
    cmd = [
        'java', '-jar', trimmomatic_path, 'PE',
        '-threads', '4',
        input_r1, input_r2,
        output_r1_paired, output_r1_unpaired,
        output_r2_paired, output_r2_unpaired,
        adapter_arg,
        'LEADING:3',
        'TRAILING:3',
        'SLIDINGWINDOW:4:15',
        'MINLEN:36'
    ]
    
    # Run Trimmomatic
    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully processed {base_name}")
        return [output_r1_paired, output_r2_paired]  # Return paths of paired output files
    except subprocess.CalledProcessError as e:
        print(f"Error processing {base_name}: {e}")
        return []

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run Trimmomatic and FastQC on multiple FASTQ files using a schema.')
    parser.add_argument('-j', '--jar', 
                       required=True,
                       help='Path to trimmomatic.jar')
    parser.add_argument('-s', '--schema',
                       default='schema_trim.tsv',
                       help='Path to schema TSV file (default: schema_trim.tsv)')
    parser.add_argument('-a', '--adapters',
                       default='adapters.fa',
                       help='Path to adapters FASTA file (default: adapters.fa)')
    parser.add_argument('--skip-fastqc',
                       action='store_true',
                       help='Skip running FastQC on trimmed files')
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Verify the jar file exists
    if not os.path.exists(args.jar):
        print(f"Error: Cannot find Trimmomatic at {args.jar}")
        return
    
    # Verify the schema file exists
    if not os.path.exists(args.schema):
        print(f"Error: Cannot find schema file at {args.schema}")
        return
        
    # Verify the adapters file exists
    if not os.path.exists(args.adapters):
        print(f"Error: Cannot find adapters file at {args.adapters}")
        return
    
    # Create QC output directory
    qc_dir = 'QC'
    os.makedirs(qc_dir, exist_ok=True)
    
    # Read the schema file
    schema = pd.read_csv(args.schema, sep='\t', header=None)
    schema.columns = ['fastq_file', 'adapter_type']
    
    # Store all trimmed file paths
    trimmed_files = []
    
    # Process each file in the schema
    for idx, row in schema.iterrows():
        print(f"Processing {row['fastq_file']} with {row['adapter_type']}")
        output_files = run_trimmomatic(row['fastq_file'], row['adapter_type'], args.adapters, args.jar)
        trimmed_files.extend(output_files)
    
    # Run FastQC on trimmed files
    if not args.skip_fastqc:
        print("\nRunning FastQC on trimmed files...")
        for trimmed_file in trimmed_files:
            run_fastqc(trimmed_file, qc_dir)

if __name__ == "__main__":
    main()
