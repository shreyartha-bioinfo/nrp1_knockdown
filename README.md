# RNA-seq Analysis of NRP1 and VEGFB in Myocardial Infarction

This repository contains scripts to analyze the correlation and co-expression network of genes associated with `NRP1` and `VEGFB` in bulk RNA-seq studies comparing Myocardial Infarction (MI) vs Control.

## Identified Studies

3 Human Studies:
- GSE249812 (Human PBMCs)
- GSE103182 (Human Whole blood)
- GSE154294 (Human Cardiomyocytes)

3 Mouse Studies:
- GSE267256 (Mouse Cardiac fibroblasts)
- GSE285626 (Mouse Cardiomyocytes)
- GSE304427 (Mouse Heart tissues)

## Scripts

- `find_studies.py`: Identifies the chosen GSE IDs.
- `fix_data.py`: Downloads the supplementary raw count/FPKM matrices from GEO.
- `parse_data.py`: Parses the raw data, mapping Ensembl IDs to gene symbols where necessary, and computes mean expressions per gene.
- `correlate_real.py`: Calculates Pearson correlations across the parsed datasets to find genes commonly moving with both NRP1 and VEGFB.
- `visualize_real.py`: Uses `networkx`, `seaborn`, and `matplotlib` to generate co-expression networks, heatmaps, and scatter plots.

## Outputs

- `correlation_results.csv`: Table of genes and their correlation scores.
- `network_diagram.png`: Network plot of NRP1, VEGFB, and top correlated genes.
- `correlation_heatmap.png`: Heatmap of pairwise correlations.
- `scatter_NRP1_VEGFB.png`: Scatter plot showing the common movement of NRP1 and VEGFB.
