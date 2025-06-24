import pandas as pd

# Load the recombination map
df = pd.read_csv("genetic_map_hg38_withX.txt", delim_whitespace=True)  # assuming header: chr, position, rate, genetic_map

# Rename columns for clarity (optional)
df.columns = ['chr', 'pos', 'rate_cM_Mb', 'genetic_map_cM']

# Threshold
threshold = 20

# Boolean for high recombination
df['high'] = df['rate_cM_Mb'] > threshold

# Group adjacent high values into clusters
df['cluster'] = (df['high'] != df['high'].shift()).cumsum()
df['cluster'] = df.apply(lambda row: row['cluster'] if row['high'] else None, axis=1)

# Select max in each cluster
local_maxima = (
    df[df['high']]
    .groupby('cluster', dropna=True)
    .apply(lambda group: group.loc[group['rate_cM_Mb'].idxmax()])
    .reset_index(drop=True)
    .drop(columns = ['high','cluster'])
)

# Output
local_maxima.to_csv("filtered_genetic_map_hg38_withX.txt",sep='\t',index=False)
