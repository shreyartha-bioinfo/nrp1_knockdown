import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import glob
from sklearn.decomposition import PCA
from scipy.stats import zscore
import networkx.algorithms.community as nx_comm

os.makedirs('results/figures', exist_ok=True)
files = glob.glob('data/processed/*.csv')

fig_count = 0

for f in files:
    gse = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, index_col=0)

    if 'NRP1' not in df.index or 'VEGFB' not in df.index:
        continue

    nrp1_corr = df.T.corrwith(df.loc['NRP1'])
    vegfb_corr = df.T.corrwith(df.loc['VEGFB'])

    combined = (nrp1_corr + vegfb_corr) / 2
    combined = combined.sort_values(ascending=False).dropna()

    local_top_genes = combined.head(30).index.tolist()
    target_genes = ['NRP1', 'VEGFB'] + [g for g in local_top_genes if g not in ['NRP1', 'VEGFB']][:25]
    top_df = df.loc[target_genes]

    variance = top_df.var(axis=1)
    meaningful_genes = variance[variance > 0.001].index.tolist()
    if 'NRP1' not in meaningful_genes or 'VEGFB' not in meaningful_genes:
        continue

    top_df = df.loc[meaningful_genes]
    pairwise_corr = top_df.T.corr()

    # 1. Advanced Network diagram
    G = nx.Graph()
    for g in meaningful_genes:
        G.add_node(g)

    threshold = 0.5
    for i in range(len(meaningful_genes)):
        for j in range(i+1, len(meaningful_genes)):
            c = pairwise_corr.iloc[i, j]
            if abs(c) > threshold and not np.isnan(c):
                G.add_edge(meaningful_genes[i], meaningful_genes[j], weight=c)

    if len(G.edges()) > 0:
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)

        # Centrality for node sizes
        degree_dict = dict(G.degree(G.nodes()))
        node_sizes = [v * 300 + 500 for v in degree_dict.values()]

        # Community detection for colors
        communities = nx_comm.greedy_modularity_communities(G)
        color_map = {}
        colors = ['lightblue', 'lightgreen', 'orange', 'pink', 'lightgray', 'yellow']
        for i, comm in enumerate(communities):
            color = colors[i % len(colors)]
            for node in comm:
                color_map[node] = color

        node_colors = [color_map.get(node, 'lightgray') for node in G.nodes()]

        # Highlight NRP1 and VEGFB
        for i, node in enumerate(G.nodes()):
            if node in ['NRP1', 'VEGFB']:
                node_colors[i] = 'red'

        edges = G.edges(data=True)
        weights = [abs(e[2]['weight']) * 3 for e in edges]
        nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=node_sizes, font_size=10, width=weights, edge_color='gray')
        plt.title(f'Co-expression Network ({gse})\nNode size by degree, color by community (NRP1/VEGFB in red)')
        plt.savefig(f'results/figures/network_{gse}.png')
        plt.close()
        fig_count += 1

    # 2. Clustermap (Z-scored)
    if len(meaningful_genes) > 2:
        try:
            z_df = top_df.apply(zscore, axis=1, result_type='broadcast').dropna(axis=1)
            if not z_df.empty:
                g_clust = sns.clustermap(z_df, cmap='viridis', figsize=(10, 10), method='ward', metric='euclidean', z_score=None)
                g_clust.fig.suptitle(f'Expression Clustermap (Z-scored) - {gse}')
                g_clust.savefig(f'results/figures/heatmap_{gse}.png')
                plt.close()
                fig_count += 1
        except Exception as e:
            print(f"Clustermap failed for {gse}: {e}")

    # 3. PCA
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

    # 4. Jointplot for NRP1 vs VEGFB (More informative than just scatter)
    try:
        g_joint = sns.jointplot(x=df.loc['NRP1'], y=df.loc['VEGFB'], kind='reg', color='darkblue')
        g_joint.fig.suptitle(f'NRP1 vs VEGFB Expression ({gse})', y=1.02)
        g_joint.savefig(f'results/figures/scatter_{gse}_NRP1_VEGFB.png')
        plt.close()
        fig_count += 1
    except:
        pass

print(f"Generated a total of {fig_count} diverse, meaningful results/figures using local correlations.")
