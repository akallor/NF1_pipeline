import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Get the top 20 pathogenic mutations

var = pd.read_csv("variants_with_metadata.tsv",sep = "\t",header = 0)

top_mutations = var[var.super_class == 'Pathogenic'].groupby('avsnp150')['Condition'].nunique().reset_index().sort_values(by='Condition', ascending=False).head(20)
top_mutation_ids = top_mutations['avsnp150'].tolist()

# Filter data to only include top mutations
filtered_data = var[(var.avsnp150.isin(top_mutation_ids)) & (var.super_class == 'Pathogenic')]

# Create a presence matrix (1 if mutation exists in condition, 0 otherwise)
presence_matrix = pd.crosstab(filtered_data['avsnp150'], filtered_data['Condition'])

# Sort the matrix by recurrence (number of conditions per mutation)
presence_matrix['total'] = presence_matrix.sum(axis=1)
presence_matrix = presence_matrix.sort_values(by='total', ascending=False)
presence_matrix = presence_matrix.drop(columns=['total'])

# Convert to binary (in case there are multiple counts)
binary_matrix = (presence_matrix > 0).astype(int)

# Create the clustered heatmap with aspect ratio control for square cells
plt.figure(figsize=(12, 12))

# Calculate aspect ratio to ensure square cells
num_rows, num_cols = binary_matrix.shape
aspect_ratio = num_cols / num_rows

# Create the heatmap with viridis colormap and no clustering (since we manually sorted)
sns.heatmap(binary_matrix,
            cmap='viridis',
            linewidths=0.5,
            cbar_kws={'label': 'Present'},
            square=True)  # Set square=True for square cells

plt.title('Heatmap of Pathogenic Mutations Across Conditions', fontsize=16)
plt.xlabel('Condition', fontsize=12)
plt.ylabel('Mutation ID (avsnp150)', fontsize=12)
plt.tight_layout()
plt.show()

# Save the figure
plt.savefig('pathogenic_mutations_heatmap_improved.pdf', bbox_inches='tight')
