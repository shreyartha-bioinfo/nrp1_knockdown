import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import glob
from sklearn.decomposition import PCA

os.makedirs('results/figures', exist_ok=True)
files = glob.glob('data/processed/*.csv')

fig_count = 0

for f in files:
    gse = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, index_col=0)

    if 'NRP1' not in df.index or 'VEGFB' not in df.index:
        continue

    # INSTEAD OF USING A GLOBAL TOP WHICH IS DOMINATED BY ONE SPECIES/DATASET
    # WE CALCULATE THE TOP CORRELATED GENES IN *THIS* SPECIFIC DATASET

    nrp1_corr = df.T.corrwith(df.loc['NRP1'])
    vegfb_corr = df.T.corrwith(df.loc['VEGFB'])

    # Combined score
    combined = (nrp1_corr + vegfb_corr) / 2
    combined = combined.sort_values(ascending=False).dropna()

    # Take top 30
    local_top_genes = combined.head(30).index.tolist()

    target_genes = ['NRP1', 'VEGFB'] + [g for g in local_top_genes if g not in ['NRP1', 'VEGFB']][:25]
    top_df = df.loc[target_genes]

    variance = top_df.var(axis=1)
    meaningful_genes = variance[variance > 0.001].index.tolist()
    if 'NRP1' not in meaningful_genes or 'VEGFB' not in meaningful_genes:
        continue

    top_df = df.loc[meaningful_genes]
    pairwise_corr = top_df.T.corr()

    # 1. Network diagram
    G = nx.Graph()
    for g in meaningful_genes:
        G.add_node(g)

    # Dynamic threshold to ensure edges exist
    threshold = 0.5
    for i in range(len(meaningful_genes)):
        for j in range(i+1, len(meaningful_genes)):
            c = pairwise_corr.iloc[i, j]
            if abs(c) > threshold and not np.isnan(c):
                G.add_edge(meaningful_genes[i], meaningful_genes[j], weight=c)

    if len(G.edges()) > 0:
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
        node_colors = ['red' if node in ['NRP1', 'VEGFB'] else 'lightblue' for node in G.nodes()]
        edges = G.edges(data=True)
        weights = [abs(e[2]['weight']) * 2 for e in edges]
        nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=1500, font_size=9, width=weights, edge_color='gray')
        plt.title(f'Network for {gse} (Threshold > 0.5)')
        plt.savefig(f'results/figures/network_{gse}.png')
        plt.close()
        fig_count += 1

    # 2. Heatmap
    if len(meaningful_genes) > 2:
        plt.figure(figsize=(12, 10))
        sns.heatmap(pairwise_corr, cmap='coolwarm', vmin=-1, vmax=1)
        plt.title(f'Correlation Heatmap - {gse}')
        plt.tight_layout()
        plt.savefig(f'results/figures/heatmap_{gse}.png')
        plt.close()
        fig_count += 1

    # 3. PCA on the subset of top genes
    pca = PCA(n_components=2)
    try:
        comps = pca.fit_transform(top_df.T)
        plt.figure(figsize=(6, 5))
        plt.scatter(comps[:, 0], comps[:, 1], alpha=0.7, color='purple')
        plt.title(f'PCA of top correlated genes ({gse})')
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} var)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} var)')
        plt.savefig(f'results/figures/pca_{gse}.png')
        plt.close()
        fig_count += 1
    except:
        pass

    # 4. Boxplots of NRP1 & VEGFB expression distributions
    plt.figure(figsize=(6, 5))
    plot_data = pd.DataFrame({'Expression': pd.concat([df.loc['NRP1'], df.loc['VEGFB']]),
                              'Gene': ['NRP1']*df.shape[1] + ['VEGFB']*df.shape[1]})
    sns.violinplot(x='Gene', y='Expression', data=plot_data, inner="point")
    plt.title(f'Expression Distribution in {gse}')
    plt.savefig(f'results/figures/violin_{gse}.png')
    plt.close()
    fig_count += 1

    # 5. Scatter Plot NRP1 vs VEGFB
    plt.figure(figsize=(6, 5))
    sns.regplot(x=df.loc['NRP1'], y=df.loc['VEGFB'])
    plt.title(f'NRP1 vs VEGFB in {gse}')
    plt.savefig(f'results/figures/scatter_{gse}_NRP1_VEGFB.png')
    plt.close()
    fig_count += 1

    # 6. Scatter plots for high-variance top genes
    high_var_genes = [g for g in meaningful_genes if g not in ['NRP1', 'VEGFB']][:4]
    for g in high_var_genes:
        plt.figure(figsize=(5, 4))
        sns.regplot(x=df.loc['NRP1'], y=df.loc[g], color='green')
        plt.title(f'{gse}: NRP1 vs {g}')
        plt.savefig(f'results/figures/scatter_{gse}_NRP1_vs_{g}.png')
        plt.close()
        fig_count += 1

        plt.figure(figsize=(5, 4))
        sns.regplot(x=df.loc['VEGFB'], y=df.loc[g], color='orange')
        plt.title(f'{gse}: VEGFB vs {g}')
        plt.savefig(f'results/figures/scatter_{gse}_VEGFB_vs_{g}.png')
        plt.close()
        fig_count += 1

print(f"Generated a total of {fig_count} diverse, meaningful results/figures using local correlations.")
