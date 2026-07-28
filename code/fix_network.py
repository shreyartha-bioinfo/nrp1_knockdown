import pandas as pd
corrs = pd.read_csv('correlation_results.csv', index_col=0)
print(corrs.head())
