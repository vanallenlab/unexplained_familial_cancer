import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load the data from the TSV files
hwe_df = pd.read_csv('/Users/noah/Desktop/ufc_repository/figures/supplemental/step_7_hwe/tp_prob_hwe_output_snp.tsv', sep='\t')
sensitivity_df = pd.read_csv('/Users/noah/Desktop/ufc_repository/figures/supplemental/step_6_sensitivity/sensitivity_output_snp.tsv', sep='\t')
denovo_df = pd.read_csv('/Users/noah/Desktop/ufc_repository/figures/supplemental/step_5_count_variants/denovo_snp.counts.tsv', sep='\t')
total_df = pd.read_csv('/Users/noah/Desktop/ufc_repository/figures/supplemental/step_5_count_variants/snp.counts.tsv', sep='\t')

# Rename column if it exists
denovo_df = denovo_df.rename(columns={"TP_PROB": "tp_prob","Q2":'median_DeNovo_SNPs_variants'})
total_df = total_df.rename(columns={"TP_PROB": "tp_prob","Q2":'median_SNPs_variants'})

# Merge all dataframes on 'tp_prob'
merged_df = pd.merge(hwe_df, sensitivity_df, on='tp_prob')
merged_df = pd.merge(merged_df, denovo_df, on='tp_prob')
merged_df = pd.merge(merged_df, total_df, on='tp_prob')

# Normalize the 'total_snps' column by dividing by the max value
merged_df['median_SNPs_norm'] = abs(merged_df['median_SNPs_variants'] - 3800000) / 3800000
merged_df['median_DeNovo_SNPs_norm'] = abs(merged_df['median_DeNovo_SNPs_variants'] - 65) / 65
merged_df['sensitivity_norm'] = merged_df['sensitivity'] / merged_df['sensitivity'].max()
merged_df['median_p_hwe_1_percent_euro_norm'] = merged_df['median_p_hwe_1_percent_euro'] / merged_df['median_p_hwe_1_percent_euro'].max()

# Define columns to optimize
max_cols = ['median_SNPs_norm', 'sensitivity_norm']  # Columns to maximize
min_cols = ['median_DeNovo_SNPs', 'median_p_hwe_1_percent_euro']  # Columns to minimize

#merged_df.head()
pd.set_option('display.max_rows', None)
# Calculate the Euclidean distance for each row
merged_df['euclidean_distance'] = np.sqrt(
    (1 * (1 - merged_df['median_SNPs_norm'])) ** 2 +  # Maximizing
    (1 * (1 - merged_df['sensitivity_norm'])) ** 2 +  # Maximizing
    (1 * merged_df['median_DeNovo_SNPs_norm']) ** 2 + # Minimizing
    (1 * merged_df['median_p_hwe_1_percent_euro_norm']) ** 2  # Minimizing
)
merged_df.to_csv("snps_tmp.tsv",sep='\t',index=False)
# Assuming your DataFrame is named df and has columns 'tp_prob' and 'euclidean_distance'
plt.plot(merged_df['tp_prob'], merged_df['euclidean_distance'], marker='o', linestyle='-', color='b')

# Adding labels and title
plt.xlabel('TP Probability')
plt.ylabel('Euclidean Distance')
plt.title('TP Probability vs Euclidean Distance (SNPs)')
plt.savefig('SNPs_TP_Threshold_Euclidean.png')