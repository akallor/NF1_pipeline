import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.signal import find_peaks


# Load data (already using the real data)
var = pd.read_csv("variants_with_metadata.tsv", sep="\t")

# Filter to ensure positions are positive
var = var[var['Start'] > 0].copy()

# Sort by chromosome and position - update to handle chrM properly
var['chr_num'] = var['Chr'].str.replace('chr', '').replace('M', '26').replace('X', '24').replace('Y', '25').astype(str)
var['chr_num'] = pd.to_numeric(var['chr_num'], errors='coerce')
var = var.sort_values(['chr_num', 'Start'])

# Filter to focus on pathogenic variants
path_variants = var[var['super_class'] == 'Pathogenic'].copy()

# Function to create the density plot for a single chromosome
def plot_chromosome_density(chromosome, ax=None, highlight_peaks=True):
    chrom_data = path_variants[path_variants['Chr'] == chromosome]

    if len(chrom_data) < 5:  # Skip if too few variants
        return None, []

    # Get positions
    positions = chrom_data['Start'].values

    # Calculate KDE (density)
    positions_range = np.linspace(min(positions), max(positions), 1000)
    kde = stats.gaussian_kde(positions)
    density = kde(positions_range)

    # Plot the density curve
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 3))

    # Plot base density
    ax.plot(positions_range, density, color='steelblue', alpha=0.8)
    ax.fill_between(positions_range, density, color='steelblue', alpha=0.3)

    # Find and highlight peaks if requested
    peak_positions = []
    if highlight_peaks and len(density) > 3:
        # Find peaks with prominence to avoid minor fluctuations
        peaks, _ = find_peaks(density, prominence=0.1*max(density))

        if len(peaks) > 0:
            # Get the top 3 peaks or fewer if less exist
            top_peak_indices = peaks[np.argsort(density[peaks])[-min(3, len(peaks)):]]

            # Highlight the peaks
            ax.plot(positions_range[top_peak_indices], density[top_peak_indices], 'ro')

            # Add annotations for peak positions (in genomic coordinates)
            for idx in top_peak_indices:
                peak_pos = positions_range[idx]
                ax.annotate(f"{int(peak_pos):,}",
                           (peak_pos, density[idx]),
                           xytext=(5, 10),
                           textcoords='offset points',
                           arrowprops=dict(arrowstyle='->', color='red'),
                           color='darkred',
                           fontsize=8)
                peak_positions.append(int(peak_pos))

    # Set labels
    ax.set_title(f"Pathogenic Variant Density - {chromosome}")
    ax.set_xlabel("Position")
    ax.set_ylabel("Density")

    # Format y-axis to be more readable
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))

    # Remove top and right spines for cleaner look
    sns.despine(ax=ax)

    return ax, peak_positions

# Create a multi-panel visualization for all chromosomes
def plot_all_chromosomes(chromosomes_to_plot=None):
    if chromosomes_to_plot is None:
        # Get list of chromosomes with at least 5 pathogenic variants
        chrom_counts = path_variants['Chr'].value_counts()
        chromosomes_to_plot = chrom_counts[chrom_counts >= 5].index.tolist()

    # Sort chromosomes naturally - update to handle chrM properly
    chromosomes_to_plot = sorted(chromosomes_to_plot,
                                key=lambda x: int(x.replace('chr', '').replace('M', '26').replace('X', '24').replace('Y', '25')))

    n_chroms = len(chromosomes_to_plot)

    if n_chroms == 0:
        print("No chromosomes with sufficient pathogenic variants to plot")
        return None, {}

    # Calculate grid dimensions: try to make it roughly square
    n_cols = min(3, n_chroms)
    n_rows = (n_chroms + n_cols - 1) // n_cols

    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 3*n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])  # Ensure axes is array-like
    axes = axes.flatten()

    # Store peak information
    peak_info = {}

    # Plot each chromosome
    for i, chrom in enumerate(chromosomes_to_plot):
        if i < len(axes):
            _, peaks = plot_chromosome_density(chrom, ax=axes[i])
            if peaks:
                peak_info[chrom] = peaks

    # Hide unused subplots
    for i in range(n_chroms, len(axes)):
        axes[i].axis('off')

    plt.tight_layout()
    return fig, peak_info

# Execute the visualization
fig, peak_info = plot_all_chromosomes()

# Print information about the highest peaks for each chromosome
if peak_info:
    print("Highest density peaks of pathogenic variants by chromosome:")
    for chrom, peaks in peak_info.items():
        if peaks:
            positions_str = ", ".join([f"{pos:,}" for pos in peaks])
            print(f"{chrom}: {positions_str}")

    plt.savefig('variant_density_distribution_1.pdf')
    plt.show()

# Create a genome-wide view with normalization by chromosome length
def plot_genome_overview():
    # Create a dataframe with chromosome lengths
    chrom_lengths = {
        'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
        'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
        'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
        'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
        'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
        'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrM': 16569, 'chrY': 59373566
    }

    # Map variants to their relative position in genome
    var_plot = path_variants.copy()
    var_plot['chrom_number'] = var_plot['Chr'].str.replace('chr', '').replace('M', '26').replace('X', '24').replace('Y', '25')
    var_plot['chrom_number'] = pd.to_numeric(var_plot['chrom_number'], errors='coerce')

    # Sort chromosomes
    var_plot = var_plot.sort_values('chrom_number')

    # Group by chromosome and count
    chrom_counts = var_plot['Chr'].value_counts()

    # Create DataFrame for plotting
    plot_data = pd.DataFrame({
        'Chromosome': chrom_counts.index,
        'Count': chrom_counts.values
    })

    # Add chromosome lengths and normalize counts
    plot_data['Length'] = plot_data['Chromosome'].map(chrom_lengths)
    plot_data['Normalized'] = (plot_data['Count'] / plot_data['Length']) * 1e6  # Variants per million bp

    # Sort chromosomes properly for display
    plot_data['Sorting'] = plot_data['Chromosome'].str.replace('chr', '').replace('M', '26').replace('X', '24').replace('Y', '25')
    plot_data['Sorting'] = pd.to_numeric(plot_data['Sorting'], errors='coerce')
    plot_data = plot_data.sort_values('Sorting')

    # Create figures (two subplots side by side)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

    # Plot 1: Raw counts
    sns.barplot(x='Chromosome', y='Count', data=plot_data, palette='viridis', ax=ax1)
    ax1.set_title('Raw Pathogenic Variants by Chromosome')
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Number of Pathogenic Variants')
    ax1.tick_params(axis='x', rotation=45)

    # Add text labels for counts
    for i, count in enumerate(plot_data['Count']):
        if count > 0:
            ax1.text(i, count + 0.5, int(count), ha='center')

    # Plot 2: Normalized counts
    sns.barplot(x='Chromosome', y='Normalized', data=plot_data, palette='rocket', ax=ax2)
    ax2.set_title('Pathogenic Variants per Million Base Pairs')
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Variants per Mb')
    ax2.tick_params(axis='x', rotation=45)

    # Add text labels for normalized counts
    for i, norm in enumerate(plot_data['Normalized']):
        if norm > 0:
            ax2.text(i, norm + 0.05, f"{norm:.2f}", ha='center')

    plt.tight_layout()
    plt.savefig('variants_density_distribution_2.pdf')
    plt.show()

plot_genome_overview()
