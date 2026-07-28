import pandas as pd
import glob
import os
import warnings
warnings.filterwarnings('ignore')

files = glob.glob('data/processed/*.csv')

all_corrs = []

for f in files:
    gse = os.path.basename(f).split('.')[0]
    df = pd.read_csv(f, index_col=0)

    if 'NRP1' in df.index and 'VEGFB' in df.index:
        nrp1_corr = df.T.corrwith(df.loc['NRP1'])
        vegfb_corr = df.T.corrwith(df.loc['VEGFB'])

        corr_df = pd.DataFrame({
            'NRP1_corr': nrp1_corr,
            'VEGFB_corr': vegfb_corr
        })
        corr_df['Study'] = gse
        corr_df['Gene'] = corr_df.index
        all_corrs.append(corr_df)
    else:
        print(f"Skipping {gse}, missing NRP1 or VEGFB")

combined = pd.concat(all_corrs)

avg_corr = combined.groupby('Gene')[['NRP1_corr', 'VEGFB_corr']].mean()
# To find genes commonly moving with *both*, we can look at the minimum of the two absolute correlations
# or just the average, but let's do the mean of absolute correlations to find strong movement
avg_corr['Combined_Score'] = (avg_corr['NRP1_corr'] + avg_corr['VEGFB_corr']) / 2

avg_corr = avg_corr.sort_values(by='Combined_Score', ascending=False)
avg_corr = avg_corr.dropna()
avg_corr.to_csv('results/correlation_results.csv')

print("Top 10 correlated genes with both NRP1 and VEGFB (Real Data):")
print(avg_corr.head(10))
