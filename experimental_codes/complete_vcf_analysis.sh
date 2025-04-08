#!/bin/bash
#SBATCH --job-name=vcf_analysis
#SBATCH --output=vcf_analysis_%j.log
#SBATCH --error=vcf_analysis_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G


# Print job info
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Start time: $(date)"



# Activate virtual environment if needed
# source /path/to/your/virtualenv/bin/activate

# Create a directory for output in the user's scratch space
OUTDIR=$PWD/vcf_analysis_${SLURM_JOB_ID}
mkdir -p $OUTDIR
#cd $OUTDIR

# Copy your input files to the compute node (modify path as needed)
#rsync -avz $PWD/JH-2-* .

# Create the Python script from your code
cat > vcf_analysis.py << 'EOF'
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np
from cyvcf2 import VCF  # For parsing VCF files efficiently

"""
Can be run for multiple samples and has been tested for one and more than one sample
"""

def extract_patient_id(filepath):
    """Extract patient ID from filepath."""
    # Extract the directory name (JH-2-XXX) from the filepath
    dirname = os.path.basename(os.path.dirname(filepath))
    # Verify it matches our pattern and return
    if dirname.startswith('JH-2-'):
        return dirname
    return None

def safe_float_conversion(value):
    """Safely convert a value to float, handling tuples and arrays."""
    if value is None:
        return None

    # If it's a tuple or list or array, take the first element
    if isinstance(value, (list, tuple, np.ndarray)):
        if len(value) > 0:
            return float(value[0])
        return None

    # Regular conversion
    try:
        return float(value)
    except (ValueError, TypeError):
        print(f"Warning: Could not convert {value} of type {type(value)} to float")
        return None

def process_vcf_files(base_dir='.', debug=False):
    """Process all VCF files with the pattern JH-2-*/*.unfiltered.vcf"""
    # Find all matching files
    pattern = os.path.join(base_dir, 'JH-2-*', '*.unfiltered.vcf')
    vcf_files = glob.glob(pattern)

    if not vcf_files:
        print(f"No files found matching pattern: {pattern}")
        return None

    print(f"Found {len(vcf_files)} VCF files to process")

    # Initialize data storage
    data = []

    # Track field types for debugging
    if debug:
        field_types = {
            'DP': set(),
            'TLOD': set(),
            'MBQ': set(),
            'BaseQRankSum': set(),
            'MQRankSum': set(),
            'BQ': set()
        }

    # Process each file
    for vcf_file in vcf_files:
        patient_id = extract_patient_id(vcf_file)
        if not patient_id:
            continue

        print(f"Processing {vcf_file}...")

        try:
            # Open VCF file using cyvcf2
            vcf = VCF(vcf_file)

            # Get the available INFO fields and print them if debugging
            if debug:
                print(f"Available INFO fields: {vcf.header_iter()}")

            # Extract values from each variant
            variant_count = 0
            for variant in vcf:
                variant_count += 1

                if debug and variant_count <= 5:
                    print(f"\nVariant {variant_count} INFO fields:")
                    info_fields = vcf.header_iter()
                    for field in info_fields:
                        if field['HeaderType'] == 'INFO':
                            field_id = field['ID']
                            value = variant.INFO.get(field_id)
                            print(f"  {field_id}: {value} (type: {type(value)})")

                # Get values from INFO fields
                dp = variant.INFO.get('DP')
                tlod = variant.INFO.get('TLOD')
                mbq = variant.INFO.get('MBQ')

                # Alternative base quality fields to try if MBQ is not found
                if mbq is None:
                    # Try other possible base quality fields
                    for field in ['BaseQRankSum', 'MQRankSum', 'BQ']:
                        mbq = variant.INFO.get(field)
                        if mbq is not None:
                            if debug:
                                print(f"Using {field} as base quality field")
                            break

                # Track types for debugging
                if debug:
                    for field, value in [('DP', dp), ('TLOD', tlod), ('MBQ', mbq)]:
                        if value is not None:
                            field_types[field].add(str(type(value)))

                # Convert values safely
                dp_float = safe_float_conversion(dp)
                tlod_float = safe_float_conversion(tlod)
                mbq_float = safe_float_conversion(mbq)

                # Store values if at least some metrics are present
                record = {'patient_id': patient_id}

                if dp_float is not None:
                    record['DP'] = dp_float

                if tlod_float is not None:
                    record['TLOD'] = tlod_float

                if mbq_float is not None:
                    record['BaseQuality'] = mbq_float

                # Add record if we have at least one metric
                if len(record) > 1:  # More than just patient_id
                    data.append(record)

            print(f"Processed {variant_count} variants from {vcf_file}")

        except Exception as e:
            print(f"Error processing {vcf_file}: {str(e)}")
            import traceback
            traceback.print_exc()

    # Print field type information if debugging
    if debug and field_types:
        print("\nField types found:")
        for field, types in field_types.items():
            print(f"{field}: {', '.join(types)}")

    # Convert to DataFrame
    if data:
        df = pd.DataFrame(data)
        print(f"Created DataFrame with {len(df)} rows and columns: {', '.join(df.columns)}")
        return df
    else:
        print("No data extracted from VCF files")
        return None

def plot_histograms(df):
    """Create histograms for DP, TLOD, and BaseQuality."""
    if df is None or df.empty:
        return

    # Determine which metrics are available
    available_metrics = [col for col in ['DP', 'TLOD', 'BaseQuality'] if col in df.columns]
    num_metrics = len(available_metrics)

    if num_metrics == 0:
        print("No metrics to plot")
        return

    # Set up the figure
    fig, axes = plt.subplots(1, num_metrics, figsize=(6*num_metrics, 5))
    if num_metrics == 1:
        axes = [axes]  # Make it iterable

    # Create histograms for each available metric
    for i, metric in enumerate(available_metrics):
        sns.histplot(data=df, x=metric, kde=True, ax=axes[i])
        axes[i].set_title(f'Distribution of {metric} Values')
        axes[i].set_xlabel(metric)
        axes[i].set_ylabel('Count')

    plt.tight_layout()
    plt.savefig('metrics_histograms.png', dpi=300)
    plt.close()

    # Create log-scale versions for metrics that might span orders of magnitude
    for metric in ['TLOD', 'DP']:
        if metric in df.columns:
            plt.figure(figsize=(8, 5))
            sns.histplot(data=df, x=metric, kde=True)
            plt.xscale('log')
            plt.title(f'Distribution of {metric} Values (Log Scale)')
            plt.xlabel(f'{metric} - Log Scale')
            plt.ylabel('Count')
            plt.tight_layout()
            plt.savefig(f'{metric.lower()}_histogram_log.png', dpi=300)
            plt.close()

def plot_boxplots(df):
    """Create boxplots for metrics by patient."""
    if df is None or df.empty:
        return

    # Extract numeric part from patient ID (assuming format JH-2-XXX)
    df['patient_num'] = df['patient_id'].apply(lambda x: int(re.search(r'JH-2-(\d+)', x).group(1))
                                              if re.search(r'JH-2-(\d+)', x) else 0)
    sorted_patients = df.sort_values('patient_num')['patient_id'].unique()

    # Determine which metrics are available
    available_metrics = [col for col in ['DP', 'TLOD', 'BaseQuality'] if col in df.columns]

    # Create boxplots for each metric
    for metric in available_metrics:
        plt.figure(figsize=(max(10, len(sorted_patients)/2), 6))
        sns.boxplot(data=df, x='patient_id', y=metric, order=sorted_patients)
        plt.title(f'{metric} Values by Patient')
        plt.xlabel('Patient ID')
        plt.ylabel(metric)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(f'{metric.lower()}_boxplot_by_patient.png', dpi=300)
        plt.close()

        # Log scale version for metrics that might span orders of magnitude
        if metric in ['TLOD', 'DP']:
            plt.figure(figsize=(max(10, len(sorted_patients)/2), 6))
            sns.boxplot(data=df, x='patient_id', y=metric, order=sorted_patients)
            plt.yscale('log')
            plt.title(f'{metric} Values by Patient (Log Scale)')
            plt.xlabel('Patient ID')
            plt.ylabel(f'{metric} - Log Scale')
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.savefig(f'{metric.lower()}_boxplot_by_patient_log.png', dpi=300)
            plt.close()

def create_correlation_plots(df):
    """Create correlation plots between metrics."""
    if df is None or df.empty:
        return

    # Determine which metrics are available
    available_metrics = [col for col in ['DP', 'TLOD', 'BaseQuality'] if col in df.columns]

    # Need at least 2 metrics to create correlation
    if len(available_metrics) < 2:
        return

    # Create pairwise scatter plots
    for i in range(len(available_metrics)):
        for j in range(i+1, len(available_metrics)):
            metric1 = available_metrics[i]
            metric2 = available_metrics[j]

            plt.figure(figsize=(8, 6))
            sns.scatterplot(data=df, x=metric1, y=metric2, alpha=0.5)
            plt.title(f'Correlation between {metric1} and {metric2}')
            plt.xlabel(metric1)
            plt.ylabel(metric2)
            plt.tight_layout()
            plt.savefig(f'{metric1.lower()}_{metric2.lower()}_correlation.png', dpi=300)
            plt.close()

            # Calculate and print correlation coefficient
            corr = df[[metric1, metric2]].corr().iloc[0, 1]
            print(f"Correlation between {metric1} and {metric2}: {corr:.4f}")

def create_summary_stats(df):
    """Create summary statistics for all metrics."""
    if df is None or df.empty:
        return

    # Determine which metrics are available
    available_metrics = [col for col in ['DP', 'TLOD', 'BaseQuality'] if col in df.columns]

    if not available_metrics:
        return

    # Overall summary
    overall_summary = df[available_metrics].describe()
    print("\nOverall Summary Statistics:")
    print(overall_summary)

    # By patient summary
    patient_summary = df.groupby('patient_id')[available_metrics].describe()
    print("\nSummary Statistics by Patient:")
    print(patient_summary)

    # Save to CSV
    overall_summary.to_csv('overall_metrics_stats.csv')
    patient_summary.to_csv('patient_metrics_stats.csv')

    return overall_summary, patient_summary

def create_heatmap(df):
    """Create heatmap of median metric values by patient."""
    if df is None or df.empty:
        return

    # Determine which metrics are available
    available_metrics = [col for col in ['DP', 'TLOD', 'BaseQuality'] if col in df.columns]

    if not available_metrics:
        return

    # Calculate median values for each patient and metric
    median_values = df.groupby('patient_id')[available_metrics].median().reset_index()

    # Extract patient number for sorting
    median_values['patient_num'] = median_values['patient_id'].apply(
        lambda x: int(re.search(r'JH-2-(\d+)', x).group(1)) if re.search(r'JH-2-(\d+)', x) else 0
    )

    # Sort by patient number
    median_values = median_values.sort_values('patient_num')

    # Prepare data for heatmap
    heatmap_data = median_values.set_index('patient_id')[available_metrics]

    # Create separate heatmaps for each metric since they might have different scales
    for metric in available_metrics:
        plt.figure(figsize=(max(10, len(heatmap_data)/2), 6))
        sns.heatmap(heatmap_data[[metric]].T, annot=True, cmap="YlGnBu", fmt=".1f",
                   linewidths=.5, cbar_kws={"label": f"Median {metric} Value"})
        plt.title(f'Median {metric} Values by Patient')
        plt.tight_layout()
        plt.savefig(f'{metric.lower()}_heatmap.png', dpi=300)
        plt.close()

    # Create combined Z-score heatmap to compare metrics on the same scale
    if len(available_metrics) > 1:
        # Calculate Z-scores for each metric
        z_scores = pd.DataFrame(index=heatmap_data.index)
        for metric in available_metrics:
            z_scores[metric] = (heatmap_data[metric] - heatmap_data[metric].mean()) / heatmap_data[metric].std()

        plt.figure(figsize=(max(10, len(z_scores)/2), 6))
        sns.heatmap(z_scores.T, annot=True, cmap="coolwarm", fmt=".2f",
                   linewidths=.5, cbar_kws={"label": "Z-Score"}, center=0)
        plt.title('Standardized Metric Values by Patient (Z-Scores)')
        plt.tight_layout()
        plt.savefig('metrics_z_score_heatmap.png', dpi=300)
        plt.close()

def main(base_dir='.'):
    """Main function to process files and create plots."""
    # Process all VCF files with debug mode enabled
    df = process_vcf_files(base_dir, debug=True)

    if df is not None and not df.empty:
        print(f"Successfully extracted data: {len(df)} records from {df['patient_id'].nunique()} patients")

        # Check which metrics were found
        available_metrics = [col for col in ['DP', 'TLOD', 'BaseQuality'] if col in df.columns]
        print(f"Available metrics: {', '.join(available_metrics)}")

        # Create plots
        plot_histograms(df)
        plot_boxplots(df)
        create_correlation_plots(df)
        create_heatmap(df)

        # Generate summary statistics
        create_summary_stats(df)

        # Save processed data
        df.to_csv('vcf_metrics_data.csv', index=False)

        print("Analysis complete. Output files saved to current directory.")
    else:
        print("No data to analyze.")

if __name__ == "__main__":
    # Get input directory from environment variable
    import os
    base_directory = os.environ.get("INPUT_DIR", ".")
    print(f"Using input directory: {base_directory}")
    main(base_directory)
EOF

# Install required Python packages
#pip install --user pandas matplotlib seaborn numpy cyvcf2

# Set the input directory where your VCF files are located
# Edit this path to match your data location
export INPUT_DIR="$SLURM_SUBMIT_DIR"

# Run the Python script
echo "Starting VCF analysis..."
python vcf_analysis.py

# Copy results back to a permanent storage location
mkdir -p $PWD/vcf_analysis_results_${SLURM_JOB_ID}
cp -r *.png *.csv $PWD/vcf_analysis_results_${SLURM_JOB_ID}/

echo "Analysis complete. Results copied to $HOME/vcf_analysis_results_${SLURM_JOB_ID}/"
echo "End time: $(date)"
