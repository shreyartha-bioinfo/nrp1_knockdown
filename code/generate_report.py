import os
import glob

# Ensure figures exist
figures = sorted(glob.glob('results/figures/*.png'))

with open('results/report.md', 'w') as f:
    f.write('# Comprehensive Analysis Report: NRP1 and VEGFB in Myocardial Infarction\n\n')
    f.write('This report compiles all diverse, meaningful visualizations generated from analyzing bulk RNA-seq data across 3 human and 3 mouse GEO studies of Myocardial Infarction.\n\n')

    datasets = set()
    for fig in figures:
        basename = os.path.basename(fig).replace('.png', '')
        parts = basename.split('_')
        for p in parts:
            if 'GSE' in p:
                datasets.add(p)

    for gse in sorted(datasets):
        f.write(f'## Dataset Analysis: {gse}\n\n')
        # PCA
        for fig in figures:
            if gse in fig and 'pca' in fig:
                f.write(f'### PCA\n![{os.path.basename(fig)}](figures/{os.path.basename(fig)})\n\n')
        # Clustered Heatmap
        for fig in figures:
            if gse in fig and 'heatmap' in fig:
                f.write(f'### Z-Scored Expression Clustermap\n![{os.path.basename(fig)}](figures/{os.path.basename(fig)})\n\n')
        # Advanced Network
        for fig in figures:
            if gse in fig and 'network' in fig:
                f.write(f'### Co-expression Network (Community & Centrality)\n![{os.path.basename(fig)}](figures/{os.path.basename(fig)})\n\n')
        # Jointplot
        for fig in figures:
            if gse in fig and 'scatter' in fig and 'NRP1_VEGFB' in fig:
                f.write(f'### NRP1 vs VEGFB Joint Distribution\n![{os.path.basename(fig)}](figures/{os.path.basename(fig)})\n\n')
        f.write('\n\n')

print("Report generated at results/report.md")
