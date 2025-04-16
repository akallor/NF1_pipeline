import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


var = pd.read_csv("variants_with_metadata.tsv",sep = "\t",header = 0)

# Get the top 20 pathogenic mutations
top_mutations = var[var.super_class == 'Pathogenic'].groupby('avsnp150')['Condition'].nunique().reset_index().sort_values(by='Condition', ascending=False).head(20)
top_mutation_ids = top_mutations['avsnp150'].tolist()

# Filter data to only include top mutations
filtered_data = var[(var.avsnp150.isin(top_mutation_ids)) & (var.super_class == 'Pathogenic')]

# Create a presence matrix (1 if mutation exists in condition, 0 otherwise)
presence_matrix = pd.crosstab(filtered_data['avsnp150'], filtered_data['Condition'])
presence_matrix = presence_matrix.reindex(top_mutation_ids)

# Convert to binary (in case there are multiple counts)
binary_matrix = (presence_matrix > 0).astype(int)

# Create the clustered heatmap
plt.figure(figsize=(12, 10))
sns.clustermap(binary_matrix,
               cmap='YlOrRd',
               linewidths=0.5,
               figsize=(12, 10),
               row_cluster=True,
               col_cluster=True,
               yticklabels=True,
               xticklabels=True,
               cbar_kws={'label': 'Present'})

plt.suptitle('Clustered Heatmap of Pathogenic Mutations Across Conditions', y=1.02, fontsize=16)
plt.tight_layout()
plt.show()

# Save the figure
plt.savefig('pathogenic_mutations_heatmap.pdf', bbox_inches='tight')
