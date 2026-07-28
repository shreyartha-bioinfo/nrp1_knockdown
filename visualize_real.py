import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

corrs = pd.read_csv('correlation_results.csv', index_col=0)
top_genes = corrs.head(15).index.tolist()

# Let's create an adjacency matrix using average correlation across all studies for these top genes
# Since reading all data frames and computing pairwise can be heavy, let's just do it for one robust study, or average it.
# We'll use GSE285626 as it has many samples and genes.

df = pd.read_csv('data_processed/GSE285626.csv', index_col=0)

# Filter top genes that exist in this df
valid_top_genes = [g for g in top_genes + ['NRP1', 'VEGFB'] if g in df.index]
top_df = df.loc[valid_top_genes]

pairwise_corr = top_df.T.corr()

G = nx.Graph()
for g in valid_top_genes:
    G.add_node(g)

threshold = 0.5
for i in range(len(valid_top_genes)):
    for j in range(i+1, len(valid_top_genes)):
        g1 = valid_top_genes[i]
        g2 = valid_top_genes[j]
        c = pairwise_corr.loc[g1, g2]
        if abs(c) > threshold:
            G.add_edge(g1, g2, weight=c)

plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G, seed=42)

node_colors = []
for node in G.nodes():
    if node in ['NRP1', 'VEGFB']:
        node_colors.append('red')
    else:
        node_colors.append('lightblue')

edges = G.edges(data=True)
weights = [abs(e[2]['weight']) * 2 for e in edges]

nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=2000,
        font_size=10, font_weight='bold', edge_color='gray', width=weights)

plt.title('Co-expression Network (Nodes: Top correlated with NRP1 & VEGFB)')
plt.savefig('network_diagram.png')
plt.close()

plt.figure(figsize=(12, 10))
sns.heatmap(pairwise_corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1, fmt=".2f")
plt.title('Correlation Heatmap of NRP1, VEGFB and Top Related Genes')
plt.tight_layout()
plt.savefig('correlation_heatmap.png')
plt.close()

plt.figure(figsize=(8, 6))
sns.regplot(x=df.loc['NRP1'], y=df.loc['VEGFB'], scatter_kws={'alpha':0.6}, color='purple')
plt.title('Common Movement: NRP1 vs VEGFB Expression')
plt.xlabel('NRP1 Expression')
plt.ylabel('VEGFB Expression')
plt.savefig('scatter_NRP1_VEGFB.png')
plt.close()

print("Real data visualizations saved.")
