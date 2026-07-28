import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
from scipy.stats import zscore

os.makedirs('results/figures', exist_ok=True)
files = glob.glob('data/processed/*.csv')

for f in files:
    gse = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, index_col=0)

    meta_file = f'data/metadata/{gse}_labels.csv'
    if not os.path.exists(meta_file):
        continue

    meta = pd.read_csv(meta_file)
    meta = meta[meta['Condition'] != 'Unknown']

    if len(meta['Condition'].unique()) < 2:
        print(f"Skipping {gse}, not enough classes.")
        continue

    # check intersection of samples. The CSV from processed data might have different columns depending on original file formats
    # Sometimes it's #12012 instead of GSM...
    # For GSE304427 the columns are WT_Con1 not GSM... so metadata mapping from soft file fails here
    # The user asked for "where are the MI v control labels". Our metadata extractor mapped to GSMs.
    # The processed CSVs have different headers because parse_data.py didn't replace them with GSMs.
    print(f"Data columns for {gse}: {df.columns[:5]}")
    print(f"Meta samples for {gse}: {meta['Sample'].values[:5]}")
