import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load the data from the TSV files
hwe_df = pd.read_csv('hwe_snps.tsv', sep='\t')
sensitivity_df = pd.read_csv('sensitivity_snps.tsv', sep='\t')
denovo_df = pd.read_csv('denovo_snps.tsv', sep='\t')
total_df = pd.read_csv('total_snps.tsv', sep='\t')

# Merge all dataframes on 'tp_prob'
merged_df = pd.merge(hwe_df, sensitivity_df, on='tp_prob')
merged_df = pd.merge(merged_df, denovo_df, on='tp_prob')
merged_df = pd.merge(merged_df, total_df, on='tp_prob')

# Normalize the 'total_snps' column by dividing by the max value
merged_df['median_SNPs_norm'] = merged_df['median_SNPs_variants'] / merged_df['median_SNPs_variants'].max()
merged_df['median_DeNovo_SNPs_norm'] = merged_df['median_DeNovo_SNPs_variants'] / merged_df['median_DeNovo_SNPs_variants'].max()
merged_df['median_p_hwe_1_percent_euro_norm'] = merged_df['median_p_hwe_1_percent_euro'] / merged_df['median_p_hwe_1_percent_euro'].max()

# Define columns to optimize
max_cols = ['median_SNPs_norm', 'sensitivity']  # Columns to maximize
min_cols = ['median_DeNovo_SNPs_variants', 'median_p_hwe_1_percent_euro']  # Columns to minimize

#merged_df.head()
pd.set_option('display.max_rows', None)
# Calculate the Euclidean distance for each row
merged_df['euclidean_distance'] = np.sqrt(
    (1-merged_df['median_SNPs_norm']) ** 2 + (1- merged_df['sensitivity']) ** 2 + (merged_df['median_DeNovo_SNPs_norm']) ** 2 + (merged_df['median_p_hwe_1_percent_euro_norm']) ** 2
)

# Assuming your DataFrame is named df and has columns 'tp_prob' and 'euclidean_distance'
plt.plot(merged_df['tp_prob'], merged_df['euclidean_distance'], marker='o', linestyle='-', color='b')

# Adding labels and title
plt.xlabel('TP Probability')
plt.ylabel('Euclidean Distance')
plt.title('TP Probability vs Euclidean Distance (SNPs)')
plt.save_fig('SNPs_TP_Threshold_Euclidean.png')